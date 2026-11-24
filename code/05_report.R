# =============================================================================
# 05_report.R  --  Writes outputs/run_report.md
# A self-contained description of this run: data, choices, settings, figures and
# an overview of the CSV. Run after 04_export.R.
# =============================================================================
source("code/00_config.R")

L   <- character()
add <- function(...) L <<- c(L, ...)

# ---- helpers ---------------------------------------------------------------------------
md_table <- function(df, digits = 3) {
  fmt <- function(x) {
    out <- if (is.numeric(x)) format(round(x, digits), trim = TRUE, scientific = FALSE)
    else as.character(x)
    out[is.na(x)] <- "NA"
    out
  }
  cells <- unname(lapply(df, fmt))
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    if (nrow(df)) paste0("| ", do.call(paste, c(cells, sep = " | ")), " |"))
}

legend_tbl <- tryCatch(terra::cats(terra::rast(CLC_PATH))[[1]], error = function(e) NULL)
code_label <- function(codes) {
  nm <- if (!is.null(legend_tbl) && all(c("LEVEL4", "Level4Eng") %in% names(legend_tbl)))
    legend_tbl$Level4Eng[match(codes, legend_tbl$LEVEL4)] else rep("", length(codes))
  paste0("`", codes, "` ", nm)
}
pkg_ver <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                                error = function(e) "not installed")
yn <- function(x) if (isTRUE(x)) "yes" else "no"
pct <- function(x) paste0(100 * x, " %")

# ---- read results ------------------------------------------------------------------------
tbl  <- read.csv(file.path(OUT_DIR, "grid_200m_all.csv"), stringsAsFactors = FALSE)
kept <- tbl[!tbl$exclude_water, ]
g    <- make_grid()
ee   <- as.vector(g$ext)                       # xmin, xmax, ymin, ymax
tif  <- terra::rast(file.path(OUT_DIR, "grid_200m.tif"))
n_all <- nrow(tbl); n_excl <- sum(tbl$exclude_water); n_kept <- nrow(kept)
n_nocost <- sum(is.na(kept$total_time_min))
fi   <- file.info(CLC_PATH)
osm_f <- file.path(CACHE_DIR, "osm_roads_route_area.gpkg")
osm_date <- if (file.exists(osm_f)) format(file.info(osm_f)$mtime, "%Y-%m-%d") else "unknown"

# =============================================================================
add("# Run report: forest, edge, water and travel-cost grid", "",
    paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), " by `05_report.R`. ",
           "Everything below is written from the settings in `00_config.R` and the files in `outputs/`."),
    "")

add("## 1. What this produces", "",
    paste0("A grid of ", N_CELLS, " x ", N_CELLS, " cells of ", CELL_M, " m (",
           N_CELLS * CELL_M / 1000, " x ", N_CELLS * CELL_M / 1000,
           " km) around Ylistö campus, Jyväskylä. For each cell: forest cover, forest edge density, ",
           "forest cover in a ", BUF_RADIUS_M / 1000, " km circle around the cell centre, ",
           "water share, and the cost of reaching it from the campus (car + walking). ",
           "Cells with more than ", pct(WATER_EXCLUDE_THRESHOLD), " open water are flagged and set to NA in the GeoTIFF."),
    "")

add("## 2. Study area and grid", "",
    paste0("- Centre: campus at ", CENTRE_LAT, " N, ", CENTRE_LON, " E (WGS84)"),
    paste0("- Grid centre snapped to the ", CELL_M, " m lattice (EPSG:3067): E ",
           g$snapped_centre[1], ", N ", g$snapped_centre[2]),
    paste0("- Extent (EPSG:3067): xmin ", ee[1], ", xmax ", ee[2], ", ymin ", ee[3], ", ymax ", ee[4]),
    paste0("- Cell size ", CELL_M, " m = ", (CELL_M / CLC_RES_M)^2, " Corine pixels (", CELL_M / CLC_RES_M,
           " x ", CELL_M / CLC_RES_M, "), so cells nest exactly"),
    paste0("- Cell centroids are given in EPSG:3067 (`x_3067`, `y_3067`) and WGS84 (`lon`, `lat`)"),
    paste0("- Extra land cover read for edges: ", EDGE_BUFFER_M, " m around the grid"),
    paste0("- Roads / hiking area: grid + ", ROUTE_BUFFER_M / 1000, " km on every side"),
    "")

add("## 3. Data used", "")
add(md_table(data.frame(
  Dataset = c("Corine Land Cover 2024", "OpenStreetMap roads"),
  Source = c("Finnish Environment Institute (SYKE)",
             "OpenStreetMap contributors, Geofabrik Finland extract via `osmextract`"),
  Detail = c(paste0("20 m raster, EPSG:3067, file `", basename(CLC_PATH), "` (", round(fi$size / 1e6),
                    " MB). LEVEL4 codes."),
             paste0("Vector lines; clipped copy cached on ", osm_date, ". Routing with `dodgr`.")),
  UsedFor = c("Forest, water, marsh, edge density, walking friction, bridges",
              "Driving network, driving time and distance, drivable bridges")
)), "")
add("Licences: check the SYKE terms for the 2024 release; OpenStreetMap data are ODbL.", "",
    paste0("Considered but not used: SYKE 2 m land cover (no tree-height classes), ",
           "Luke multi-source forest inventory 16 m, a DEM (no terrain slope in this version)."),
    "")

add("## 4. Definitions and choices", "", "### 4.1 Forest", "",
    "A 20 m Corine pixel is forest if its LEVEL4 code is one of:", "",
    paste0("- ", code_label(FOREST_CODES)), "",
    paste0("Everything else, including transitional woodland/shrub (3241-3246), agro-forestry (2441) ",
           "and urban green (1411), is non-forest. The Corine legend gives canopy-cover thresholds but ",
           "no tree height, so a 5 m height rule cannot be checked directly."),
    "",
    "- `forest_prop` = forest pixels / all pixels of the cell",
    "- `forest_prop_land` = forest pixels / non-water pixels of the cell", "")

add("### 4.2 Edge density", "",
    paste0("Edge = boundary between a forest pixel and a non-forest pixel (4-neighbour adjacency: ",
           "left/right/up/down, no diagonals). Each differing pair counts ", CLC_RES_M, " m. ",
           "Half of it is given to each of the two pixels, so an edge along a cell border is shared ",
           "between the two cells. Edge length per cell is divided by the cell area (",
           CELL_M^2 / 1e4, " ha) to give m/ha."),
    "",
    "- Forest vs ALL non-forest counts as edge, including forest-water shorelines (`edge_density_m_ha`)",
    "- `edge_density_nowater_m_ha` ignores forest-water boundaries",
    "- Counting along pixel sides overestimates diagonal boundaries and misses features narrower than one 20 m pixel; use it to rank cells, not as an absolute value to compare with other studies",
    "")

add(paste0("### 4.2b Forest cover in a ", BUF_RADIUS_M / 1000, " km buffer"), "",
    paste0("`", BUF_COL, "` is the share of forest pixels in a circle of radius ", BUF_RADIUS_M,
           " m centred on the centre of the ", CELL_M, " m cell. It is counted directly on the ",
           CLC_RES_M, " m Corine pixels, using the same forest codes as `forest_prop`."),
    "",
    paste0("- A pixel is inside the circle if its centre is within ", BUF_RADIUS_M, " m of the cell centre (",
           format(round(4 * sum(floor((sqrt(4 * (BUF_RADIUS_M / CLC_RES_M)^2 -
                                            (2 * seq_len(BUF_RADIUS_M / CLC_RES_M) - 1)^2) + 1) / 2))),
                  big.mark = ","), " pixels); pixels are counted whole, none are split"),
    "- Denominator = all pixels in the circle: water and pixels without data count as non-forest (like `forest_prop`), and the focal cell is part of the circle",
    "- `buf_valid_prop` = share of the circle with Corine data (1 unless the circle reaches outside the raster)",
    "- Computed before the water exclusion: cells flagged for water still contribute to the circles of their neighbours",
    "")

add("### 4.3 Water and exclusion", "",
    "Open water = Corine codes:", "", paste0("- ", code_label(WATER_CODES)), "",
    paste0("Aquatic marsh is not counted as water. A cell is excluded (`exclude_water = TRUE`) if ",
           "`water_prop` > ", WATER_EXCLUDE_THRESHOLD, ". Excluded cells are kept in `grid_200m_all.csv`, ",
           "removed from `grid_200m_filtered.csv`, NA in `grid_200m.tif` and grey in the figures. ",
           "The cost surface is computed before exclusion because lakes shape the routes."),
    "", "Marsh / bog codes (`marsh_prop`, slow walking):", "", paste0("- ", code_label(MARSH_CODES)), "")

add("### 4.4 Driving", "",
    paste0("- Network: OpenStreetMap, car profile of `dodgr` (fastest route), area = grid + ",
           ROUTE_BUFFER_M / 1000, " km"),
    paste0("- Origin: the road vertex nearest to the campus"),
    paste0("- Road types kept: ", paste(HIGHWAY_KEEP, collapse = ", ")),
    paste0("- Forest roads (`highway=track` and forestry/agricultural access) allowed: ", yn(ALLOW_FOREST_ROADS),
           if (ALLOW_FOREST_ROADS) paste0(" (grades ", paste(FOREST_TRACK_GRADES, collapse = ", "),
                                          "; driven at ", FOREST_ROAD_SPEED_KMH, " km/h)") else ""),
    paste0("- Private roads allowed: ", yn(ALLOW_PRIVATE_ROADS), "; access = no is always excluded; ",
           "the most specific tag wins (motorcar > motor_vehicle > access)"),
    "- Speeds: OSM numeric `maxspeed` where present, otherwise dodgr's defaults for the road class",
    "- Distance is measured along the fastest route, which can be a little longer than the shortest route",
    "- Not modelled: winter conditions, gates, ferries, traffic, junction delays beyond dodgr's defaults", "")

add("### 4.5 Hiking", "",
    paste0("- Working grid ", HIKE_RES_M, " m, 8-neighbour moves; a diagonal move needs both orthogonal neighbours passable"),
    paste0("- Walking speed on normal ground: ", WALK_SPEED_KMH, " km/h"),
    paste0("- Open water: impassable. A ", HIKE_RES_M, " m cell is passable if at least ",
           pct(PASSABLE_MIN_FRACTION), " of its 20 m pixels are passable (narrow rivers can vanish)"),
    paste0("- Marsh / bog / wet ground: ",
           if (MARSH_IMPASSABLE) "impassable" else paste0("walked at ", MARSH_SPEED_KMH, " km/h")),
    "- Other slow ground (walking time multiplier):",
    paste0("  - ", code_label(FRICTION_OVERRIDES$code), ": x", FRICTION_OVERRIDES$friction),
    paste0("- Drivable bridges over water treated as walkable ground: ", yn(BRIDGES_PASSABLE),
           " (only roads that are in the car routing graph)"),
    "- No terrain slope in this version", "")

add("### 4.6 How a cell's travel cost is found", "",
    paste0("Every road vertex is a possible parking/access point with its drive time and distance. ",
           "For each cell the pipeline picks the access point that minimises:"),
    "",
    paste0("`", OBJ_PER_DRIVE_MIN, " x drive minutes + ", OBJ_PER_DRIVE_KM, " x drive km + ",
           OBJ_PER_HIKE_MIN, " x hike minutes`"),
    "",
    paste0("The default (1, 0, 1) means fastest total trip. Hiking runs from that access point to the ",
           "hiking cell under the cell centroid (or the nearest reachable one inside the cell; ",
           "`cost_offset_m` gives the offset). Path distance is the length of the least-time walking route, ",
           "not a straight line. All drive/hike values are one way."),
    "")

add("### 4.7 Fuel and money columns", "",
    paste0("- Consumption ", FUEL_L_PER_100KM, " L/100 km, price ", FUEL_PRICE_EUR_PER_L, " EUR/L, value of time ",
           VALUE_OF_TIME_EUR_H, " EUR/h, round trip: ", yn(ROUND_TRIP)),
    "- `fuel_l`, `fuel_cost_eur`, `time_cost_eur`, `total_cost_eur` are computed from `drive_dist_km` and the times (doubled if round trip)",
    "")

# ---- results -----------------------------------------------------------------------------------
add("## 5. Results", "",
    paste0("- Cells: ", n_all, " in total; ", n_excl, " excluded for water > ", pct(WATER_EXCLUDE_THRESHOLD),
           "; ", n_kept, " kept"),
    paste0("- Kept cells without a travel-cost value: ", n_nocost), "")
if ("cost_status" %in% names(kept)) {
  st <- as.data.frame(table(kept$cost_status, useNA = "ifany"))
  names(st) <- c("cost_status (kept cells)", "cells")
  add(md_table(st), "")
}
if (file.exists(file.path(OUT_DIR, "grid_landcover_overview.png")))
  add("### Land cover", "", "![Forest, edge density and water](grid_landcover_overview.png)", "")
if (file.exists(file.path(OUT_DIR, "grid_cost_surface.png")))
  add("### Cost surface", "", "![Driving and hiking cost](grid_cost_surface.png)", "")
add("Grey cells = water > limit. Colour scale: ", paste0(PALETTE_NAME, if (PALETTE_REVERSED) " (reversed)" else ""), ".", "")

# ---- CSV overview ----------------------------------------------------------------------------------
add("## 6. Output files", "",
    paste0("- `grid_200m_all.csv`: ", n_all, " rows x ", ncol(tbl), " columns (all cells, `exclude_water` flag)"),
    paste0("- `grid_200m_filtered.csv`: ", n_kept, " rows (water <= ", pct(WATER_EXCLUDE_THRESHOLD), ")"),
    paste0("- `grid_200m.tif`: ", terra::ncol(tif), " x ", terra::nrow(tif), " cells, ", terra::nlyr(tif),
           " layers (", paste(names(tif), collapse = ", "), "); excluded cells are NA"),
    "- `grid_landcover_overview.png`, `grid_cost_surface.png`, this report", "")

dict <- c(
  cell_id = "Integer cell number 1..N, row-major from the top-left cell",
  cell_uid = "Stable ID: easting/northing (EPSG:3067, m) of the cell's south-west corner",
  row = "Grid row (1 = top)", col = "Grid column (1 = left)",
  x_3067 = "Cell centroid easting, EPSG:3067 (m)", y_3067 = "Cell centroid northing, EPSG:3067 (m)",
  lon = "Centroid longitude, WGS84", lat = "Centroid latitude, WGS84",
  forest_prop = "Forest share of the whole cell (0-1)",
  forest_prop_land = "Forest share of the non-water area",
  edge_density_m_ha = "Forest / non-forest boundary, m per ha (shorelines count)",
  edge_density_nowater_m_ha = "Same, ignoring forest-water boundaries",
  water_prop = "Open-water share of the cell (0-1)", marsh_prop = "Marsh / bog share of the cell",
  valid_prop = "Share of the cell with valid land-cover data",
  exclude_water = "TRUE if water_prop exceeds the limit",
  drive_dist_km = "Driving distance, campus to the access point used (one way)",
  drive_time_min = "Driving time to the access point (one way)",
  hike_time_min = "Walking time from the access point to the cell (one way)",
  hike_dist_km = "Walking distance along the walking route (one way)",
  total_time_min = "Drive + hike time", total_dist_km = "Drive + hike distance",
  access_objective = "Value of the quantity minimised when choosing the access point",
  cost_offset_m = "Distance from the cell centroid to the sampled hiking cell",
  cost_status = "ok / no_passable_ground / unreachable_from_road",
  fuel_l = "Fuel use (doubled if round trip)", fuel_cost_eur = "Fuel cost, EUR",
  time_cost_eur = "Time cost, EUR (0 unless a value of time is set)",
  total_cost_eur = "Fuel + time cost, EUR")
dict[BUF_COL] <- paste0("Forest share of the ", BUF_RADIUS_M, " m circle around the cell centre (0-1; water counts as non-forest)")
dict["buf_valid_prop"] <- "Share of that circle with valid land-cover data"
cols_present <- names(tbl)
add("### Columns", "")
add(md_table(data.frame(column = cols_present,
                        meaning = ifelse(cols_present %in% names(dict), dict[cols_present], ""))), "")

skip <- c("cell_id", "row", "col", "x_3067", "y_3067", "lon", "lat")
num_cols <- setdiff(names(kept)[vapply(kept, is.numeric, logical(1))], skip)
sm <- function(f, z) suppressWarnings(f(kept[[z]], na.rm = TRUE))
summ <- data.frame(column = num_cols,
                   min    = sapply(num_cols, sm, f = min),
                   median = sapply(num_cols, sm, f = median),
                   mean   = sapply(num_cols, sm, f = mean),
                   max    = sapply(num_cols, sm, f = max),
                   n_NA   = sapply(num_cols, function(z) sum(is.na(kept[[z]]))),
                   row.names = NULL)
add("### Summary of kept cells", "", md_table(summ), "")

head_cols <- intersect(c("cell_id", "cell_uid", "lon", "lat", "forest_prop", "edge_density_m_ha", BUF_COL,
                         "water_prop", "drive_dist_km", "drive_time_min", "hike_time_min",
                         "hike_dist_km", "total_time_min"), names(kept))
add("### First rows (kept cells)", "", md_table(utils::head(kept[, head_cols], 6)), "")

# ---- caveats and reproducibility ---------------------------------------------------------------------
add("## 7. Known limits", "",
    "- Corine is a 20 m product built from several sources: narrow roads, power lines and small clearings are under-represented, so absolute edge densities are low. Check SYKE's documentation for the reference years of its inputs.",
    "- OSM speed limits are incomplete; defaults fill the gaps. Times assume normal summer driving and walking, no slope.",
    "- Walking speeds, marsh speed and the passable-cell rule are assumptions, not measurements.",
    "- Results have not been checked against field data or a second land-cover product.", "")

add("## 8. Reproduce", "",
    "Run `source(\"run_all.R\")` from the pipeline folder after editing `00_config.R`. Settings used:", "",
    "```r")
cfg <- c("CELL_M", "N_CELLS", "EDGE_BUFFER_M", "ROUTE_BUFFER_M", "BUF_RADIUS_M", "FOREST_CODES", "WATER_CODES",
         "MARSH_CODES", "WATER_EXCLUDE_THRESHOLD", "WALK_SPEED_KMH", "HIKE_RES_M", "MARSH_IMPASSABLE",
         "MARSH_SPEED_KMH", "BRIDGES_PASSABLE", "BRIDGE_FRICTION", "PASSABLE_MIN_FRACTION",
         "HIGHWAY_KEEP", "ALLOW_PRIVATE_ROADS", "ALLOW_FOREST_ROADS", "FOREST_TRACK_GRADES",
         "FOREST_ROAD_SPEED_KMH", "OBJ_PER_DRIVE_MIN", "OBJ_PER_DRIVE_KM", "OBJ_PER_HIKE_MIN",
         "FUEL_L_PER_100KM", "FUEL_PRICE_EUR_PER_L", "VALUE_OF_TIME_EUR_H", "ROUND_TRIP")
add(vapply(cfg, function(nm) paste0(nm, " <- ", paste(deparse(get(nm)), collapse = " ")), ""))
add("```", "",
    paste0("R ", R.version.string, "; terra ", pkg_ver("terra"), ", sf ", pkg_ver("sf"),
           ", dodgr ", pkg_ver("dodgr"), ", igraph ", pkg_ver("igraph"),
           ", osmextract ", pkg_ver("osmextract")))

out <- file.path(OUT_DIR, "run_report.md")
writeLines(enc2utf8(L), out, useBytes = TRUE)
message("Report written to ", normalizePath(out))
