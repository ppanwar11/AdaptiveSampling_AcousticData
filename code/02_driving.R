# =============================================================================
# 02_driving.R
# Driving time and distance from Ylistö campus to every road vertex within the
# study grid + 10 km buffer (OpenStreetMap roads, dodgr routing, no Java needed).
# Output: cache/drive_vertices.rds
# =============================================================================
source("code/00_config.R")
library(osmextract)
library(dodgr)
t0 <- Sys.time()

g <- make_grid()
ext_route <- buffer_ext(g$ext, ROUTE_BUFFER_M)

# ---- 1. Routing area as a lon/lat polygon ---------------------------------------
bb <- sf::st_bbox(c(xmin = terra::xmin(ext_route), ymin = terra::ymin(ext_route),
                    xmax = terra::xmax(ext_route), ymax = terra::ymax(ext_route)),
                  crs = sf::st_crs(3067))
area_4326 <- sf::st_transform(sf::st_as_sfc(bb), 4326)

# ---- 2. OSM roads (downloaded once, then cached) --------------------------------
# NOTE: osmextract downloads the Geofabrik extract of all Finland the first time
# (several hundred MB) and caches it; only the clipped result is kept per run.
roads_file <- file.path(CACHE_DIR, "osm_roads_route_area.gpkg")
roads <- NULL
if (file.exists(roads_file)) roads <- sf::st_read(roads_file, quiet = TRUE)
# an older cache without the tracktype tag is rebuilt (needed for forest roads)
if (is.null(roads) || !("tracktype" %in% names(roads))) {
  message("Downloading/reading OSM roads (first run is slow) ...")
  roads <- osmextract::oe_get(
    OSM_PLACE, layer = "lines",
    boundary = area_4326, boundary_type = "clipsrc",
    extra_tags = c("maxspeed", "access", "motor_vehicle", "motorcar",
                   "oneway", "junction", "surface", "tracktype"),
    query = "SELECT * FROM 'lines' WHERE highway IS NOT NULL",
    quiet = FALSE)
  sf::st_write(roads, roads_file, delete_dsn = TRUE, quiet = TRUE,
               layer = "osm_roads_route_area")
  roads <- sf::st_read(roads_file, quiet = TRUE)
}
message("  ", nrow(roads), " OSM highway segments in the routing area")

# ---- 3. Keep roads a car can use --------------------------------------------------
keep_types <- HIGHWAY_KEEP
if (ALLOW_FOREST_ROADS) keep_types <- c(keep_types, "track")
roads <- roads[roads$highway %in% keep_types, ]

# Forest tracks: only the accepted grades (when allowed at all)
if (ALLOW_FOREST_ROADS && "tracktype" %in% names(roads)) {
  is_trk <- roads$highway == "track"
  roads  <- roads[!is_trk | (roads$tracktype %in% FOREST_TRACK_GRADES), ]
}

# Access rules. In OSM the most specific tag wins: motorcar > motor_vehicle > access.
eff <- if ("access" %in% names(roads)) roads$access else rep(NA_character_, nrow(roads))
if ("motor_vehicle" %in% names(roads)) eff <- ifelse(!is.na(roads$motor_vehicle), roads$motor_vehicle, eff)
if ("motorcar"      %in% names(roads)) eff <- ifelse(!is.na(roads$motorcar),      roads$motorcar,      eff)
blocked <- "no"
if (!ALLOW_PRIVATE_ROADS) blocked <- c(blocked, "private")
if (!ALLOW_FOREST_ROADS)  blocked <- c(blocked, "forestry", "agricultural")
roads <- roads[!(eff %in% blocked), ]

# dodgr reads numeric maxspeed tags; anything else falls back to its defaults
if (!"maxspeed" %in% names(roads)) roads$maxspeed <- NA_character_
roads$maxspeed <- ifelse(grepl("^[0-9]+$", roads$maxspeed), roads$maxspeed, NA_character_)

# Tracks are not necessarily in dodgr's car profile (edges of unlisted road types are
# dropped). Treat them as an "unclassified" road with a fixed speed instead.
is_trk <- roads$highway == "track"
if (any(is_trk)) {
  roads$highway[is_trk]  <- "unclassified"
  roads$maxspeed[is_trk] <- as.character(FOREST_ROAD_SPEED_KMH)
}
roads <- suppressWarnings(sf::st_cast(roads, "LINESTRING"))
message("  ", nrow(roads), " segments kept for car routing (",
        sum(is_trk), " are forest tracks; ALLOW_FOREST_ROADS = ", ALLOW_FOREST_ROADS, ")")

# ---- 4. Routing graph -------------------------------------------------------------
graph <- dodgr::weight_streetnet(roads, wt_profile = "motorcar")
if ("component" %in% names(graph)) graph <- graph[graph$component == 1, ]  # largest network
if ("highway" %in% names(graph)) {
  message("  road types left in the routing graph (edges):")
  print(table(graph$highway))   # a type missing here was dropped by dodgr's car profile
}

# Drivable roads = only the ways that actually ended up in the routing graph (right road
# type, car access, connected to the main network). Step 03 uses these lines to find
# drivable bridges over water; footbridges, paths and closed roads never count.
drv <- roads
if ("way_id" %in% names(graph) && "osm_id" %in% names(drv)) {
  drv <- drv[as.character(drv$osm_id) %in% unique(as.character(graph$way_id)), ]
}
if (nrow(drv) == 0) { drv <- roads; warning("Could not match graph ways to roads; using all kept roads for bridges.") }
sf::st_write(sf::st_transform(drv[, "highway"], 3067),
             file.path(CACHE_DIR, "roads_drivable.gpkg"), delete_dsn = TRUE, quiet = TRUE)
message("  ", nrow(drv), " drivable road segments saved for bridge detection")

verts <- dodgr::dodgr_vertices(graph)
message("  graph: ", nrow(graph), " edges, ", nrow(verts), " vertices")

# Origin = nearest road vertex to the campus
origin_ll <- cbind(CENTRE_LON, CENTRE_LAT)
ix        <- dodgr::match_points_to_verts(verts, origin_ll)
origin_id <- verts$id[ix]
snap_m <- as.numeric(sf::st_distance(
  sf::st_transform(sf::st_sfc(sf::st_point(c(CENTRE_LON, CENTRE_LAT)), crs = 4326), 3067),
  sf::st_transform(sf::st_sfc(sf::st_point(c(verts$x[ix], verts$y[ix])), crs = 4326), 3067)))
message("  campus snapped to nearest road vertex: ", round(snap_m), " m away")
if (snap_m > 1000) warning("Origin is > 1 km from the nearest road vertex; check the road filter.")

# ---- 5. One-to-all fastest routes -------------------------------------------------
message("Routing from the campus to all vertices ...")
tm <- dodgr::dodgr_times(graph, from = origin_id, to = verts$id)                    # seconds
ds <- dodgr::dodgr_dists(graph, from = origin_id, to = verts$id, shortest = FALSE)  # metres along fastest route
ids <- as.character(verts$id)
verts$drive_time_s <- as.numeric(tm[1, ids])
verts$drive_dist_m <- as.numeric(ds[1, ids])
verts <- verts[is.finite(verts$drive_time_s), ]

# Vertex coordinates in EPSG:3067 for the hiking step
pts <- sf::st_transform(sf::st_as_sf(verts, coords = c("x", "y"), crs = 4326), 3067)
xy  <- sf::st_coordinates(pts)
out <- data.frame(id = verts$id, x3067 = xy[, 1], y3067 = xy[, 2],
                  drive_time_s = verts$drive_time_s, drive_dist_m = verts$drive_dist_m)
saveRDS(out, file.path(CACHE_DIR, "drive_vertices.rds"))

message("Reachable road vertices: ", nrow(out),
        "; max drive time ", round(max(out$drive_time_s) / 60, 1), " min, ",
        "max distance ", round(max(out$drive_dist_m) / 1000, 1), " km")
message("Done (", round(difftime(Sys.time(), t0, units = "mins"), 1), " min)")
