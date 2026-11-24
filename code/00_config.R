# =============================================================================
# 00_config.R  --  Settings and helpers shared by all scripts
# Adaptive-sampling grid: 50 x 50 km of 200 m cells around Ylistö campus,
# Jyväskylä, Finland.  Land cover: Corine Land Cover 2024 (20 m, EPSG:3067).
#
# One-time install:
#   install.packages(c("terra", "sf", "igraph", "dodgr", "osmextract"))
#
# Run everything from the folder that holds these scripts:
#   setwd("path/to/pipeline"); source("run_all.R")
# =============================================================================
setwd("/Users/popanwar/JYU/Projects/Adaptive_Sampling")

suppressPackageStartupMessages({
  library(terra)
  library(sf)
})

# ---- Paths (EDIT THIS according local project) -------------------------------------------------------
CLC_PATH  <- "data/Clc2024_FI20m/Clc2024_FI20m.tif"
OUT_DIR   <- "outputs"
CACHE_DIR <- "cache"
dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Study grid ---------------------------------------------------------------
CENTRE_LON <- 25.74144217494367      # Ylistö campus
CENTRE_LAT <- 62.22928094770784
CRS_FI     <- "EPSG:3067"           # ETRS-TM35FIN
CELL_M     <- 200                    # analysis cell size (m)
N_CELLS    <- 250                    # 250 x 200 m = 50 km per side
CLC_RES_M  <- 20                     # Corine pixel size (m)
EDGE_BUFFER_M  <- 1000               # extra land cover around grid for edge calc
ROUTE_BUFFER_M <- 10000              # extra area around grid for roads / hiking

# ---- Forest cover in a circular buffer around each cell ---------------------------
# Circle centred on the CENTRE of each 200 m cell. A 20 m Corine pixel belongs to the
# circle if its own centre is within BUF_RADIUS_M of the cell centre (about 7,850
# pixels for 1000 m). Share = forest pixels / ALL pixels of the circle, so water and
# no-data count as non-forest (same convention as forest_prop). The focal cell is part
# of the circle. Needs: BUF_RADIUS_M a multiple of CLC_RES_M, CELL_M / CLC_RES_M even,
# and EDGE_BUFFER_M + CELL_M / 2 >= BUF_RADIUS_M (checked in 01).
# The column name follows the radius: 1000 -> "forest_buf1km", 500 -> "forest_buf0.5km".
BUF_RADIUS_M <- 1000
BUF_COL      <- paste0("forest_buf", BUF_RADIUS_M / 1000, "km")

# ---- Corine 2024 classes (LEVEL4 codes, from your legend) -----------------------
# Forest = all "forest" classes (broad-leaved, coniferous, mixed; mineral soil,
# peatland, rocky).  NOT forest: transitional woodland/shrub (32xx, canopy
# cover < 30 %), agro-forestry (2441), urban green (1411).
FOREST_CODES <- c(3111, 3112, 3121, 3122, 3123, 3131, 3132, 3133)
# Open water (rivers, lakes, sea).  Aquatic marsh (4112, 4212) is NOT counted
# as water here; it is counted as marsh.
WATER_CODES  <- c(5111, 5121, 5231)
MARSH_CODES  <- c(4111, 4112, 4121, 4211, 4212)   # marshes + open peatbogs
WATER_EXCLUDE_THRESHOLD <- 0.20      # cells with water_prop > this are flagged

# ---- Hiking cost -----------------------------------------------------------------
WALK_SPEED_KMH   <- 4                # straight-line walking speed on friction 1
HIKE_RES_M       <- 40               # working resolution; must divide 200 and be a
                                     # multiple of 20 (20, 40 or 100). 20 is most
                                     # detailed but needs ~4x the memory of 40.
MARSH_IMPASSABLE <- FALSE            # TRUE = marsh/bog cannot be crossed at all
MARSH_SPEED_KMH  <- 1                # walking speed in marsh / bog / wet ground
                                     # (normal ground = WALK_SPEED_KMH)
# DRIVABLE bridges: where a road that is in the car routing network (step 02) crosses
# open water (Corine "water" pixels), that crossing is treated as solid ground. Cars
# already drive over it in step 02; this makes the bridge usable in step 03 as well, so
# its vertices can be access points and you can walk from it. Footbridges, paths and
# roads the car filter removed are NOT counted.
BRIDGES_PASSABLE <- TRUE
BRIDGE_FRICTION  <- 1                # 1 = normal walking speed on the bridge
PASSABLE_MIN_FRACTION <- 0.5         # a hiking cell is passable if >= this share
                                     # of its 20 m pixels is passable
# Other slow ground (friction multiplier on walking time). Edit freely.
FRICTION_OVERRIDES <- data.frame(
  code     = c(3112, 3122, 3132, 3243, 4122),   # peatland forests/woodland, peat production
  friction = c(1.5,  1.5,  1.5,  1.5,  1.5)
)

# ---- Driving ---------------------------------------------------------------------
OSM_PLACE <- "Finland"
HIGHWAY_KEEP <- c("motorway", "motorway_link", "trunk", "trunk_link",
                  "primary", "primary_link", "secondary", "secondary_link",
                  "tertiary", "tertiary_link", "unclassified", "residential",
                  "living_street", "service")        # forest tracks: see ALLOW_FOREST_ROADS
ALLOW_PRIVATE_ROADS <- FALSE         # FALSE = drop roads whose access is "private"
                                     # (access "no" is always dropped)

# ---- Forest roads: QUICK TOGGLE ----------------------------------------------------------
# FALSE: forest roads are NOT usable. Ways mapped highway=track are left out, as are
#        roads with access / motor_vehicle / motorcar = forestry or agricultural.
# TRUE : tracks of the grades below, and forestry/agricultural-access roads, are driven
#        at FOREST_ROAD_SPEED_KMH. Many forest roads are tagged access=private; if
#        turning this on changes little, also set ALLOW_PRIVATE_ROADS <- TRUE.
# After changing this, re-run steps 02, 03 and 04.
ALLOW_FOREST_ROADS    <- FALSE
FOREST_TRACK_GRADES   <- c("grade1", "grade2", "grade3", NA)   # NA = tracktype not mapped
FOREST_ROAD_SPEED_KMH <- 30

# ---- Which road access point is used to reach each cell? ----------------------------
# For every cell the pipeline picks the parking/access point that MINIMISES
#   OBJ_PER_DRIVE_MIN * drive_minutes + OBJ_PER_DRIVE_KM * drive_km
#   + OBJ_PER_HIKE_MIN * hike_minutes
# Units are up to you (they only set the trade-off).  Examples:
#   fastest overall (default):  1, 0, 1
#   fuel-saving first:          0.01, 1, 0.01   (kilometres dominate, minutes break ties)
#   euros:                      time_value_eur_per_min, fuel_eur_per_km, time_value_eur_per_min
# OBJ_PER_HIKE_MIN must be > 0 (otherwise hiking would be "free" and paths arbitrary).
# The driving route itself is always the fastest route (OSM car profile).
OBJ_PER_DRIVE_MIN <- 1
OBJ_PER_DRIVE_KM  <- 0
OBJ_PER_HIKE_MIN  <- 1

# ---- Fuel / money columns (computed in step 4; change and re-run 04 only) -------------
FUEL_L_PER_100KM     <- 7.0          # car consumption
FUEL_PRICE_EUR_PER_L <- 1.90         # EDIT to the current pump price
VALUE_OF_TIME_EUR_H  <- 0            # set > 0 to add a time cost (driving + hiking)
ROUND_TRIP           <- TRUE         # TRUE: fuel and time costs count there AND back

# ---- Figures and GeoTIFF ---------------------------------------------------------------
# Cells with water_prop > WATER_EXCLUDE_THRESHOLD are grey in the PNGs and NA in the TIF.
PALETTE_NAME      <- "viridis"       # any grDevices::hcl.colors() name
PALETTE_REVERSED  <- TRUE            # TRUE: low = yellow, high = dark purple (viridis);
                                     # FALSE: low = dark purple, high = yellow. Flip if
                                     # the ramp is the wrong way round for you.
EXCLUDED_COLOUR   <- "grey70"
WRITE_UNMASKED_TIF <- FALSE          # TRUE: also write a TIF where water cells are NOT NA

# =============================================================================
# Helpers
# =============================================================================

# Extent buffered by d metres on every side
buffer_ext <- function(e, d) {
  terra::ext(terra::xmin(e) - d, terra::xmax(e) + d,
             terra::ymin(e) - d, terra::ymax(e) + d)
}

# Ylistö in EPSG:3067
centre_xy <- function() {
  p <- sf::st_sfc(sf::st_point(c(CENTRE_LON, CENTRE_LAT)), crs = 4326)
  sf::st_coordinates(sf::st_transform(p, 3067))[1, c("X", "Y")]
}

# The 200 m study grid, centred on the campus and snapped to the 200 m lattice
# (so it lines up exactly with the 20 m Corine pixels).
make_grid <- function() {
  xy   <- centre_xy()
  cx   <- round(xy[["X"]] / CELL_M) * CELL_M
  cy   <- round(xy[["Y"]] / CELL_M) * CELL_M
  half <- N_CELLS * CELL_M / 2
  e    <- terra::ext(cx - half, cx + half, cy - half, cy + half)
  tmpl <- terra::rast(e, resolution = CELL_M, crs = CRS_FI)
  list(ext = e, template = tmpl, centre = xy, snapped_centre = c(cx, cy))
}

# Read a window of the Corine raster as plain numeric codes.
# Returns list(geom = empty SpatRaster with the window's geometry,
#              codes = numeric vector, row-major from the top-left)
read_clc_window <- function(window_ext) {
  clc <- terra::rast(CLC_PATH)
  r   <- terra::crop(clc, window_ext)
  v   <- terra::values(r, mat = FALSE)
  if (!is.numeric(v)) {
    stop("Corine values came back as ", class(v)[1],
         ", not numeric codes. Check terra::values() on a small crop.")
  }
  message("  window ", ncol(r), " x ", nrow(r), " pixels; code range ",
          paste(range(v, na.rm = TRUE), collapse = "-"))
  list(geom = terra::rast(r), codes = as.numeric(v))
}

# Put a vector (row-major) into a raster with the given geometry
vec_to_raster <- function(geom, x) {
  r <- terra::rast(geom)
  terra::values(r) <- x
  r
}

METRIC_NAMES <- c("forest_prop", "forest_prop_land", "edge_density_m_ha",
                  "edge_density_nowater_m_ha", "water_prop", "marsh_prop",
                  "valid_prop", BUF_COL, "buf_valid_prop")
