# =============================================================================
# 04_export.R
# Combine land-cover metrics and travel cost into the final CSV, GeoTIFF and figures.
#
# outputs/grid_200m_all.csv            all 62,500 cells, with exclude_water flag
# outputs/grid_200m_filtered.csv       only cells with water_prop <= 0.20
# outputs/grid_200m.tif                multi-band GeoTIFF; cells with water_prop > 0.20 are NA
# outputs/grid_200m_unmasked.tif       (optional, WRITE_UNMASKED_TIF) same without the NA mask
# outputs/grid_landcover_overview.png  forest / edge / water maps + histograms
# outputs/grid_cost_surface.png        driving, hiking and fuel-cost maps
# =============================================================================
source("code/00_config.R")

g <- make_grid()
n <- terra::ncell(g$template)

m <- terra::rast(file.path(CACHE_DIR, "metrics_200m.tif"))
if (terra::nlyr(m) != length(METRIC_NAMES))
  stop("cache/metrics_200m.tif has ", terra::nlyr(m), " layers but METRIC_NAMES has ",
       length(METRIC_NAMES), ". Re-run 01_landcover_metrics.R (it now adds ", BUF_COL, ").")
names(m) <- METRIC_NAMES
cost <- readRDS(file.path(CACHE_DIR, "cost_200m.rds"))
stopifnot(terra::ncell(m) == n, nrow(cost) == n)

# ---- Table -------------------------------------------------------------------------
cent <- terra::xyFromCell(g$template, seq_len(n))
rc   <- terra::rowColFromCell(g$template, seq_len(n))
ll   <- sf::st_coordinates(sf::st_transform(
          sf::st_as_sf(data.frame(x = cent[, 1], y = cent[, 2]),
                       coords = c("x", "y"), crs = 3067), 4326))
M <- terra::values(m)

tbl <- data.frame(
  cell_id = seq_len(n), row = rc[, 1], col = rc[, 2],
  x_3067 = cent[, 1], y_3067 = cent[, 2],
  lon = ll[, 1], lat = ll[, 2],
  forest_prop               = M[, "forest_prop"],
  #forest_prop_land          = M[, "forest_prop_land"],
  edge_density_m_ha         = M[, "edge_density_m_ha"],
  #edge_density_nowater_m_ha = M[, "edge_density_nowater_m_ha"],
  water_prop                = M[, "water_prop"],
  marsh_prop                = M[, "marsh_prop"],
  valid_prop                = M[, "valid_prop"]
)
tbl[[BUF_COL]]       <- M[, BUF_COL]            # forest share of the circle around the cell centre
tbl$buf_valid_prop   <- M[, "buf_valid_prop"]   # share of that circle with Corine data
tbl <- cbind(tbl, cost[, setdiff(names(cost), "cell_id")])
tbl$exclude_water <- !is.na(tbl$water_prop) & tbl$water_prop > WATER_EXCLUDE_THRESHOLD

# Fuel and money columns (driving distance/time stay one-way in drive_*; these
# include the round-trip factor if ROUND_TRIP is TRUE). Edit the prices in
# 00_config.R and re-run only this script.
trip <- if (ROUND_TRIP) 2 else 1
tbl$fuel_l         <- tbl$drive_dist_km * FUEL_L_PER_100KM / 100 * trip
tbl$fuel_cost_eur  <- tbl$fuel_l * FUEL_PRICE_EUR_PER_L
tbl$time_cost_eur  <- (tbl$drive_time_min + tbl$hike_time_min) / 60 *
                      VALUE_OF_TIME_EUR_H * trip
tbl$total_cost_eur <- tbl$fuel_cost_eur + tbl$time_cost_eur

# Column order: ids, coordinates, habitat, the travel columns you asked for, the rest
first_cols <- c("cell_id", "row", "col", "x_3067", "y_3067", "lon", "lat",
                "forest_prop", "edge_density_m_ha", BUF_COL, "water_prop", "exclude_water",
                "drive_dist_km", "drive_time_min", "hike_time_min", "hike_dist_km",
                "total_time_min", "total_dist_km")
tbl <- tbl[, c(first_cols, setdiff(names(tbl), first_cols))]

# ---- CSV -----------------------------------------------------------------------------
num <- vapply(tbl, is.numeric, logical(1))
tbl[num] <- lapply(tbl[num], function(z) round(z, 5))
write.csv(tbl, file.path(OUT_DIR, "grid_200m_all.csv"), row.names = FALSE)
write.csv(tbl[!tbl$exclude_water, ], file.path(OUT_DIR, "grid_200m_filtered.csv"),
          row.names = FALSE)

# ---- RDS (full precision, before the CSV rounding) ----------------------------------------
saveRDS(tbl, file.path(OUT_DIR, "grid_200m_table.rds"))

# ---- GeoTIFF --------------------------------------------------------------------------
value_cols <- setdiff(names(tbl), c("cell_id", "row", "col", "x_3067", "y_3067",
                                    "lon", "lat", "exclude_water", "cost_status"))
layers <- lapply(value_cols, function(nm) {
  r <- terra::rast(g$template); terra::values(r) <- tbl[[nm]]; names(r) <- nm; r
})
vals_r <- do.call(c, layers)                       # all cells, all layers

excl_r <- terra::rast(g$template)                  # 1 = excluded (water > threshold), else NA
terra::values(excl_r) <- ifelse(tbl$exclude_water, 1, NA)
masked <- terra::mask(vals_r, excl_r, maskvalues = 1)   # excluded cells -> NA in every layer

terra::writeRaster(masked, file.path(OUT_DIR, "grid_200m.tif"),
                   overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)
if (WRITE_UNMASKED_TIF) {
  terra::writeRaster(vals_r, file.path(OUT_DIR, "grid_200m_unmasked.tif"),
                     overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)
}

# ---- Figures -------------------------------------------------------------------------------
pal <- grDevices::hcl.colors(100, PALETTE_NAME, rev = PALETTE_REVERSED)
any_excl <- any(tbl$exclude_water)
thr_txt  <- paste0("water > ", 100 * WATER_EXCLUDE_THRESHOLD, " %")

# Map of one layer; excluded cells (water > threshold) are drawn grey.
draw_map <- function(r, main, grey = TRUE) {
  ok <- tryCatch({
    terra::plot(r, col = pal, main = main, mar = c(2, 2, 2.5, 4))
    TRUE
  }, error = function(e) FALSE)
  if (!ok) { plot.new(); title(main = paste(main, "(no data)")); return(invisible()) }
  if (grey && any_excl) {
    terra::plot(excl_r, col = EXCLUDED_COLOUR, legend = FALSE, add = TRUE)
    legend("bottomleft", legend = thr_txt, fill = EXCLUDED_COLOUR, bty = "n", cex = 0.8)
  }
}

draw_hist <- function(x, main, xlab, thr = NULL) {
  x <- x[is.finite(x)]
  hist(x, breaks = 30, col = pal[50], border = "white", main = main, xlab = xlab,
       ylab = "number of cells")
  abline(v = median(x), lty = 2)
  legend("topright", legend = c(paste0("n = ", length(x)),
                                paste0("median = ", signif(median(x), 3))),
         bty = "n", cex = 0.8)
  if (!is.null(thr)) abline(v = thr, col = "red", lwd = 2)
}

kept <- !tbl$exclude_water

# Figure 1: land cover
buf_txt <- paste0(BUF_RADIUS_M / 1000, " km")
png(file.path(OUT_DIR, "grid_landcover_overview.png"), width = 3200, height = 1500, res = 200)
op <- par(mfrow = c(2, 4))
draw_map(masked[["forest_prop"]],       "Forest cover (proportion of cell)")
draw_map(masked[["edge_density_m_ha"]], "Forest edge density (m/ha)")
draw_map(masked[[BUF_COL]],             paste0("Forest cover in ", buf_txt, " buffer (proportion)"))
draw_map(vals_r[["water_prop"]],        "Water (proportion of cell, all cells)", grey = FALSE)
par(mar = c(4.5, 4.5, 3, 1))
draw_hist(tbl$forest_prop[kept],       "Forest cover, kept cells",       "forest proportion")
draw_hist(tbl$edge_density_m_ha[kept], "Edge density, kept cells",       "edge density (m/ha)")
draw_hist(tbl[[BUF_COL]][kept],        paste0("Forest in ", buf_txt, " buffer, kept cells"),
          "forest proportion of circle")
draw_hist(tbl$water_prop,              "Water, all cells (red = limit)", "water proportion",
          thr = WATER_EXCLUDE_THRESHOLD)
par(op); dev.off()

# Figure 2: cost surface
# Road-network panel: line width and colour by road class (fast = thick, slow = thin).
# Uses cache/roads_drivable.gpkg written by 02_driving.R (no need to re-run 02).
ROAD_CLASSES <- list(
  list(label = "Motorway / trunk",        types = c("motorway", "motorway_link", "trunk", "trunk_link"),
       lwd = 3.0, col = "#7a0177"),
  list(label = "Primary / secondary",     types = c("primary", "primary_link", "secondary", "secondary_link"),
       lwd = 2.0, col = "#c51b8a"),
  list(label = "Tertiary / unclassified", types = c("tertiary", "tertiary_link", "unclassified"),
       lwd = 1.0, col = "#f768a1"),
  list(label = "Residential / service",   types = c("residential", "living_street", "service"),
       lwd = 0.4, col = "grey45")
)
draw_roads <- function(main = "Drivable road network (line width = road class)") {
  rf <- file.path(CACHE_DIR, "roads_drivable.gpkg")
  if (!file.exists(rf)) { plot.new(); title(main = paste(main, "(file missing)")); return(invisible()) }
  rd <- sf::st_read(rf, quiet = TRUE)
  bb <- sf::st_bbox(c(xmin = terra::xmin(g$ext), xmax = terra::xmax(g$ext),
                      ymin = terra::ymin(g$ext), ymax = terra::ymax(g$ext)), crs = sf::st_crs(3067))
  rd <- suppressWarnings(sf::st_crop(rd, bb))
  par(mar = c(2, 2, 2.5, 1))
  plot(NA, xlim = c(terra::xmin(g$ext), terra::xmax(g$ext)),
       ylim = c(terra::ymin(g$ext), terra::ymax(g$ext)), asp = 1, main = main,
       xlab = "", ylab = "", axes = TRUE)
  if (any(tbl$exclude_water)) terra::plot(excl_r, col = EXCLUDED_COLOUR, legend = FALSE, add = TRUE)
  used <- character(0)
  for (cl in rev(ROAD_CLASSES)) {                     # slow roads first, fast roads on top
    sub <- rd[rd$highway %in% cl$types, ]
    if (nrow(sub)) { plot(sf::st_geometry(sub), add = TRUE, lwd = cl$lwd, col = cl$col)
      used <- c(used, cl$label) }
  }
  points(g$centre[["X"]], g$centre[["Y"]], pch = 8, cex = 1.5, lwd = 2)
  keep <- Filter(function(cl) cl$label %in% used, ROAD_CLASSES)
  legend("bottomleft", bg = "white", cex = 0.7, title = "Road class",
         legend = c(vapply(keep, `[[`, "", "label"), "Campus (start)"),
         lwd = c(vapply(keep, `[[`, 1, "lwd"), NA), col = c(vapply(keep, `[[`, "", "col"), "black"),
         pch = c(rep(NA, length(keep)), 8))
}

# Figure 2: cost surface + road network
png(file.path(OUT_DIR, "grid_cost_surface_v2.png"), width = 3200, height = 1500, res = 200)
op <- par(mfrow = c(2, 4))
draw_roads()
draw_map(masked[["drive_time_min"]], "Driving time (min, one way)")
draw_map(masked[["drive_dist_km"]],  "Driving distance (km, one way)")
draw_map(masked[["hike_time_min"]],  "Hiking time (min, one way)")
draw_map(masked[["hike_dist_km"]],   "Hiking distance (km, one way)")
draw_map(masked[["total_time_min"]], "Total time, drive + hike (min)")
draw_map(masked[["fuel_cost_eur"]],  paste0("Fuel cost (EUR", if (ROUND_TRIP) ", round trip" else "", ")"))
par(mar = c(4.5, 4.5, 3, 1))
draw_hist(tbl$total_time_min[kept], "Total time, kept cells", "minutes (one way)")
par(op); dev.off()

# ---- Report ----------------------------------------------------------------------------------
message("Cells total: ", n, "; excluded for water > ",
        100 * WATER_EXCLUDE_THRESHOLD, " %: ", sum(tbl$exclude_water),
        "; kept: ", sum(kept))
message("Kept cells without a travel-cost value: ", sum(is.na(tbl$total_time_min[kept])))
message("Kept cells whose ", BUF_RADIUS_M, " m circle has incomplete Corine data (buf_valid_prop < 1): ",
        sum(tbl$buf_valid_prop[kept] < 1, na.rm = TRUE))
message("Correlation of the covariates (kept cells):")
print(round(cor(tbl[kept, c("forest_prop", "edge_density_m_ha", BUF_COL)],
                use = "pairwise.complete.obs"), 2))
if ("cost_status" %in% names(tbl)) print(table(tbl$cost_status[kept], useNA = "ifany"))
message("Files written to ", normalizePath(OUT_DIR))

