# ==============================================================================
# Laemanctus distribution map
# ==============================================================================
#
#
# Expected repository layout:
#
#   data/
#     Laemanctus_distribution.xlsx
#   R/
#     Laemanctus_distribution_map.R
#   output/
#
# Run this script from the repository root.
#
# Outputs:
#   output/Laemanctus_distribution_map.svg
#   output/Laemanctus_distribution_map.pdf
#   output/Laemanctus_distribution_map.png
#   output/sessionInfo.txt
#
# The first run downloads elevation data with elevatr and caches the processed
# raster in data/derived/. Subsequent runs reuse the cached GeoTIFF.
# ============================================================================== 

# ---- 1. Required packages -----------------------------------------------------

required_packages <- c(
  "sf",
  "ggplot2",
  "dplyr",
  "rnaturalearth",
  "rnaturalearthdata",
  "ggspatial",
  "elevatr",
  "terra",
  "readxl",
  "svglite"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script, for example:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(elevatr)
library(terra)
library(readxl)

# ---- 2. File paths ------------------------------------------------------------

input_file <- file.path("data", "Laemanctus_distributionSupplMat.xlsx")
output_dir <- "output"
derived_dir <- file.path("data", "derived")
elevation_file <- file.path(
  derived_dir,
  "Laemanctus_elevation_z5_WGS84.tif"
)

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nRun the script from the repository root and place the dataset in data/."
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 3. Read and validate occurrence data ------------------------------------

lizards <- readxl::read_excel(input_file)

required_columns <- c("Longitude", "Latitude", "Species")
missing_columns <- setdiff(required_columns, names(lizards))

if (length(missing_columns) > 0) {
  stop(
    "The input file is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

n_input <- nrow(lizards)

lizards_clean <- lizards %>%
  mutate(
    Longitude = as.numeric(gsub(",", ".", as.character(Longitude))),
    Latitude = as.numeric(gsub(",", ".", as.character(Latitude))),
    Species = trimws(as.character(Species))
  ) %>%
  filter(
    !is.na(Longitude),
    !is.na(Latitude),
    !is.na(Species),
    Species != ""
  )

n_removed <- n_input - nrow(lizards_clean)
if (n_removed > 0) {
  message("Removed ", n_removed, " record(s) with missing/invalid coordinates or Species.")
}

if (any(lizards_clean$Longitude < -180 | lizards_clean$Longitude > 180) ||
    any(lizards_clean$Latitude < -90 | lizards_clean$Latitude > 90)) {
  stop("At least one coordinate falls outside valid longitude/latitude ranges.")
}

expected_species <- c(
  "Laemanctus serratus",
  "Laemanctus longipes complex",
  "Laemanctus julioi",
  "Laemanctus sp."
)

unexpected_species <- setdiff(unique(lizards_clean$Species), expected_species)
missing_species <- setdiff(expected_species, unique(lizards_clean$Species))

if (length(unexpected_species) > 0) {
  stop(
    "Unexpected values in Species: ",
    paste(unexpected_species, collapse = ", ")
  )
}

if (length(missing_species) > 0) {
  warning(
    "Expected taxa absent from the dataset: ",
    paste(missing_species, collapse = ", ")
  )
}

message("Occurrence counts by taxon:")
print(table(lizards_clean$Species))

# ---- 4. Convert occurrences to sf --------------------------------------------

lizards_sf <- sf::st_as_sf(
  lizards_clean,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

# ---- 5. Map extent ------------------------------------------------------------

xmin <- -102
xmax <- -82
ymin <- 10.5
ymax <- 25

bbox <- sf::st_bbox(
  c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
  crs = sf::st_crs(4326)
)

bbox_sf <- sf::st_sf(geometry = sf::st_as_sfc(bbox))

# ---- 6. Country boundaries ----------------------------------------------------
# Keep complete country polygons. They are clipped only at plotting time by
# coord_sf(), preventing artificial border lines at the map frame.

world <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
)

# ---- 7. Elevation data --------------------------------------------------------
# The processed DEM is cached for reproducibility. Delete the cached file if you
# intentionally want to download and regenerate the elevation layer.

if (file.exists(elevation_file)) {
  message("Using cached elevation raster: ", elevation_file)
  elevation_terra <- terra::rast(elevation_file)
} else {
  message("Downloading elevation data (first run only)...")

  elevation_raw <- elevatr::get_elev_raster(
    locations = bbox_sf,
    z = 5,
    clip = "tile",
    serial = TRUE
  )

  elevation_terra <- terra::rast(elevation_raw)
  elevation_terra <- terra::project(
    elevation_terra,
    "EPSG:4326",
    method = "bilinear"
  )

  map_extent <- terra::ext(xmin, xmax, ymin, ymax)
  elevation_terra <- terra::crop(elevation_terra, map_extent)

  # Remove ocean cells while preserving complete country polygons for plotting.
  land_vector <- terra::vect(world)
  elevation_terra <- terra::mask(elevation_terra, land_vector)

  terra::writeRaster(
    elevation_terra,
    elevation_file,
    overwrite = TRUE
  )

  message("Cached processed elevation raster: ", elevation_file)
}

# Ensure the cached raster is restricted to the current map extent.
map_extent <- terra::ext(xmin, xmax, ymin, ymax)
elevation_terra <- terra::crop(elevation_terra, map_extent)
names(elevation_terra) <- "Elevation"

elevation_df <- as.data.frame(
  elevation_terra,
  xy = TRUE,
  na.rm = TRUE
)

# ---- 8. Taxon-specific layers -------------------------------------------------
# Separate layers are intentional. Laemanctus serratus is plotted last so that
# overlapping records from other taxa do not obscure it.

longipes <- lizards_sf %>%
  filter(Species == "Laemanctus longipes complex")

laemanctus_sp <- lizards_sf %>%
  filter(Species == "Laemanctus sp.")

julioi <- lizards_sf %>%
  filter(Species == "Laemanctus julioi")

serratus <- lizards_sf %>%
  filter(Species == "Laemanctus serratus")

# ---- 9. Country labels --------------------------------------------------------
# Label positions are manually adjusted for readability and to avoid dense
# occurrence clusters.

country_labels <- data.frame(
  country = c(
    "Mexico",
    "Belize",
    "Guatemala",
    "El Salvador",
    "Honduras",
    "Nicaragua",
    "Cuba"
  ),
  longitude = c(
    -97.0,
    -87.5,
    -90.8,
    -89.8,
    -86.2,
    -85.05,
    -84.0
  ),
  latitude = c(
    17.5,
    17.0,
    15.0,
    13.05,
    14.65,
    11.9,
    21.5
  )
)

# ---- 10. Elevation legend position -------------------------------------------

if (packageVersion("ggplot2") >= package_version("3.5.0")) {
  elevation_legend_theme <- theme(
    legend.position = "inside",
    legend.position.inside = c(0.22, 0.16),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = NA, color = NA),
    legend.key = element_rect(fill = NA, color = NA)
  )
} else {
  elevation_legend_theme <- theme(
    legend.position = c(0.22, 0.16),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = NA, color = NA),
    legend.key = element_rect(fill = NA, color = NA)
  )
}

# ---- 11. Build map ------------------------------------------------------------

map <- ggplot() +
  geom_tile(
    data = elevation_df,
    aes(x = x, y = y, fill = Elevation)
  ) +
  scale_fill_gradient(
    low = "gray95",
    high = "gray35",
    na.value = "transparent",
    name = "Elevation (m)",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      label.position = "bottom",
      barwidth = grid::unit(5.0, "cm"),
      barheight = grid::unit(0.30, "cm")
    )
  ) +
  geom_sf(
    data = world,
    fill = NA,
    color = "black",
    linewidth = 0.45
  ) +
  geom_sf(
    data = longipes,
    shape = 21,
    size = 3.8,
    fill = "#F2C21F",
    color = "black",
    stroke = 0.45
  ) +
  geom_sf(
    data = laemanctus_sp,
    shape = 21,
    size = 3.8,
    fill = "white",
    color = "black",
    stroke = 0.55
  ) +
  geom_sf(
    data = julioi,
    shape = 21,
    size = 3.8,
    fill = "#F44336",
    color = "black",
    stroke = 0.50
  ) +
  geom_sf(
    data = serratus,
    shape = 21,
    size = 3.9,
    fill = "#32CDB4",
    color = "black",
    stroke = 0.50
  ) +
  geom_text(
    data = country_labels,
    aes(x = longitude, y = latitude, label = country),
    size = 4.5,
    color = "black"
  ) +
  coord_sf(
    xlim = c(xmin, xmax),
    ylim = c(ymin, ymax),
    expand = FALSE,
    crs = sf::st_crs(4326),
    default_crs = sf::st_crs(4326),
    datum = sf::st_crs(4326),
    label_axes = list(
      top = "E",
      left = "N"
    ),
    clip = "on"
  ) +
  ggspatial::annotation_scale(
    location = "tr",
    width_hint = 0.17,
    unit_category = "metric",
    pad_x = grid::unit(1.55, "in"),
    pad_y = grid::unit(0.30, "in"),
    text_cex = 0.9,
    line_width = 0.8
  ) +
  ggspatial::annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = grid::unit(0.18, "in"),
    pad_y = grid::unit(0.18, "in"),
    style = ggspatial::north_arrow_orienteering(text_size = 0.8)
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.background = element_rect(
      fill = "#c6d9ec",
      color = "black",
      linewidth = 0.6
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x.top = element_text(
      size = 10,
      color = "black",
      margin = margin(b = 6)
    ),
    axis.text.x.bottom = element_blank(),
    axis.text.y.left = element_text(
      size = 10,
      color = "black",
      margin = margin(r = 6)
    ),
    axis.text.y.right = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = grid::unit(2.2, "mm"),
    axis.ticks.x.top = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x.bottom = element_blank(),
    axis.ticks.y.left = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y.right = element_blank(),
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.6
    ),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(t = 32, r = 18, b = 15, l = 26)
  ) +
  elevation_legend_theme

print(map)

# ---- 12. Export figure --------------------------------------------------------

svg_file <- file.path(output_dir, "Laemanctus_distribution_map.svg")
pdf_file <- file.path(output_dir, "Laemanctus_distribution_map.pdf")
png_file <- file.path(output_dir, "Laemanctus_distribution_map.png")

ggsave(
  filename = svg_file,
  plot = map,
  device = svglite::svglite,
  width = 12,
  height = 8,
  units = "in"
)

ggsave(
  filename = pdf_file,
  plot = map,
  device = "pdf",
  width = 12,
  height = 8,
  units = "in"
)

ggsave(
  filename = png_file,
  plot = map,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600
)

# ---- 13. Record software environment -----------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

message("Finished. Files written to: ", normalizePath(output_dir))
