# =============================================================================
# 03_hiking_cost.R
# Combined drive + hike travel cost from Ylistö campus to every 200 m cell.
#
# Method: the study area (grid + 10 km) is turned into a graph of hiking cells
# (8-neighbour moves, constant walking speed x friction; water impassable).
# A virtual "source" node is linked to every road vertex with an edge whose
# weight is that vertex's DRIVE cost (OBJ_PER_DRIVE_MIN * minutes + OBJ_PER_DRIVE_KM
# * km).  One Dijkstra run from the source then gives, for every cell, the minimum
# over all road access points of (drive cost + hiking cost from that point).
# With the default weights (1, 0, 1) this is the shortest total travel time; put
# weight on OBJ_PER_DRIVE_KM to favour low driving distance (fuel) instead.
# The access point actually used is recovered from the shortest-path tree, so
# drive and hike parts (time and distance) are reported separately.
#
# Needs: cache/drive_vertices.rds  (from 02_driving.R)
# Output: cache/cost_200m.rds, cache/cost_hike_res.tif
# =============================================================================
source("00_config.R")
library(igraph)
t0 <- Sys.time()

stopifnot(CELL_M %% HIKE_RES_M == 0, HIKE_RES_M %% CLC_RES_M == 0,
          OBJ_PER_HIKE_MIN > 0, OBJ_PER_DRIVE_MIN >= 0, OBJ_PER_DRIVE_KM >= 0,
          OBJ_PER_DRIVE_MIN + OBJ_PER_DRIVE_KM > 0)
g         <- make_grid()
ext_route <- buffer_ext(g$ext, ROUTE_BUFFER_M)

# ---- 1. Friction surface at 20 m, then aggregate to HIKE_RES_M -------------------
message("Reading Corine for the routing area ...")
win <- read_clc_window(ext_route)
v   <- win$codes

fr <- rep(1, length(v))                                   # friction 1 = normal ground
ov <- match(v, FRICTION_OVERRIDES$code)
fr[!is.na(ov)] <- FRICTION_OVERRIDES$friction[ov[!is.na(ov)]]
# Marsh / bog / wet ground: slow (e.g. 1 km/h instead of 4 km/h => friction 4)
fr[v %in% MARSH_CODES] <- if (MARSH_IMPASSABLE) NA_real_ else WALK_SPEED_KMH / MARSH_SPEED_KMH
fr[v %in% WATER_CODES] <- NA_real_                        # lakes, rivers, sea: impassable
fr[is.na(v)]           <- NA_real_                        # no data

# Bridges: drivable road lines (from step 02) that cross open-water pixels are walkable
is_br <- rep(FALSE, length(v))
road_file <- file.path(CACHE_DIR, "roads_drivable.gpkg")
if (BRIDGES_PASSABLE && file.exists(road_file)) {
  rd     <- terra::vect(road_file)
  road20 <- terra::rasterize(rd, terra::rast(win$geom), field = 1, background = 0)
  is_br  <- (terra::values(road20, mat = FALSE) == 1) & (v %in% WATER_CODES)
  is_br[is.na(is_br)] <- FALSE
  fr[is_br] <- BRIDGE_FRICTION
  message("  ", sum(is_br), " road pixels on water are treated as bridges")
}

fact   <- HIKE_RES_M / CLC_RES_M
fr20   <- vec_to_raster(win$geom, fr)
pass20 <- vec_to_raster(win$geom, as.numeric(!is.na(fr)))
br20   <- vec_to_raster(win$geom, as.numeric(is_br))
if (fact > 1) {
  pass_frac <- terra::aggregate(pass20, fact, "mean")
  fr_mean   <- terra::aggregate(fr20,   fact, "mean", na.rm = TRUE)  # mean over passable pixels
  br_h      <- terra::aggregate(br20,   fact, "max")
} else {
  pass_frac <- pass20; fr_mean <- fr20; br_h <- br20
}
# A hiking cell is passable if mostly passable, OR if a bridge crosses it (so a narrow
# bridge over a lake is not lost when the map is aggregated).
fr_h <- terra::ifel((pass_frac >= PASSABLE_MIN_FRACTION) | (br_h > 0), fr_mean, NA)

nr <- nrow(fr_h); nc <- ncol(fr_h); n <- nr * nc
Fv   <- terra::values(fr_h, mat = FALSE)                   # row-major
pass <- is.finite(Fv)
nn   <- sum(pass)
message("  hiking grid ", nc, " x ", nr, "; passable cells ", nn,
        " (", round(100 * nn / n, 1), " %)")

node_of_cell <- integer(n); node_of_cell[pass] <- seq_len(nn)
cell_of_node <- c(which(pass), NA_integer_)               # last entry = virtual source
Fm <- matrix(Fv, nr, nc, byrow = TRUE)
Im <- matrix(seq_len(n), nr, nc, byrow = TRUE)
speed <- WALK_SPEED_KMH * 1000 / 3600                      # m/s

# ---- 2. Edges between neighbouring passable cells --------------------------------
# Each undirected edge is built once, in 4 forward directions: E, SW, S, SE.
# Diagonal moves are only allowed if both orthogonal neighbours are passable,
# so a diagonal line of water cannot be "squeezed through".
message("Building hiking graph ...")
dirs <- list(c(0, 1), c(1, -1), c(1, 0), c(1, 1))
from_l <- vector("list", 4); to_l <- vector("list", 4); w_l <- vector("list", 4)
for (k in seq_along(dirs)) {
  dr <- dirs[[k]][1]; dc <- dirs[[k]][2]
  rows <- seq_len(nr - dr)
  cols <- if (dc == 1) seq_len(nc - 1) else if (dc == 0) seq_len(nc) else 2:nc
  a  <- Im[rows, cols];            b  <- Im[rows + dr, cols + dc]
  fa <- Fm[rows, cols];            fb <- Fm[rows + dr, cols + dc]
  keep <- is.finite(fa) & is.finite(fb)
  if (dr != 0 && dc != 0) {
    keep <- keep & is.finite(Fm[rows + dr, cols]) & is.finite(Fm[rows, cols + dc])
  }
  len <- HIKE_RES_M * sqrt(dr^2 + dc^2)
  from_l[[k]] <- node_of_cell[a[keep]]
  to_l[[k]]   <- node_of_cell[b[keep]]
  hike_min    <- len * (fa[keep] + fb[keep]) / 2 / speed / 60   # minutes
  w_l[[k]]    <- OBJ_PER_HIKE_MIN * hike_min                    # objective units
}

# ---- 3. Virtual source linked to road vertices (weight = drive time) -----------
rv <- readRDS(file.path(CACHE_DIR, "drive_vertices.rds"))
cell <- terra::cellFromXY(fr_h, cbind(rv$x3067, rv$y3067))
inside <- !is.na(cell)
inside[inside] <- pass[cell[inside]]
rv <- rv[inside, ]; rv$cell <- cell[inside]
rv <- rv[order(rv$drive_time_s), ]
rv <- rv[!duplicated(rv$cell), ]                           # fastest road vertex per hiking cell
src_nodes <- node_of_cell[rv$cell]
S <- nn + 1L
message("  ", length(src_nodes), " road access cells linked to the source")

from <- c(unlist(from_l), rep(S, length(src_nodes)))
to   <- c(unlist(to_l),   src_nodes)
w    <- c(unlist(w_l),
          OBJ_PER_DRIVE_MIN * rv$drive_time_s / 60 + OBJ_PER_DRIVE_KM * rv$drive_dist_m / 1000)
gr   <- igraph::make_graph(edges = as.vector(rbind(from, to)), n = S, directed = FALSE)
igraph::E(gr)$weight <- w
rm(from_l, to_l, w_l); invisible(gc())
message("  graph: ", igraph::vcount(gr), " nodes, ", igraph::ecount(gr), " edges")

# ---- 4. Shortest paths from the source ---------------------------------------------
message("Running Dijkstra ...")
dist <- igraph::distances(gr, v = S, mode = "all", weights = w)[1, ]
rm(gr); invisible(gc())

# Recover the shortest-path tree from the distances alone (no dependence on igraph
# internals): an edge a-b is a "tree edge" into b if dist[a] + weight = dist[b].
# Hiking weights are > 0, so any such parent gives a valid, loop-free tree.
parent <- integer(S)                                       # 0 = no parent
for (dirn in 1:2) {
  a  <- if (dirn == 1) from else to
  b  <- if (dirn == 1) to   else from
  da <- dist[a]; db <- dist[b]
  tight <- is.finite(da) & is.finite(db) & b != S & abs(da + w - db) <= 1e-9 * pmax(1, db)
  parent[b[tight]] <- a[tight]
}
rm(from, to, a, b, da, db, tight); invisible(gc())
pred <- parent; pred[pred == 0L] <- NA_integer_

orphans <- is.finite(dist[seq_len(nn)]) & is.na(pred[seq_len(nn)])
if (any(orphans)) stop(sum(orphans), " reachable cells have no parent in the path tree.")

is_root <- !is.na(pred) & pred == S                        # reached directly from the source
has_par <- !is.na(pred) & !is_root                         # reached from another hiking cell
up <- seq_len(S); up[has_par] <- pred[has_par]             # roots/unreached point to themselves

# Length of the last hiking step into each node
len <- numeric(S)
hp  <- which(has_par)
stopifnot(all(dist[pred[hp]] <= dist[hp] + 1e-6))          # parents are never farther away
cp  <- cell_of_node[pred[hp]]; cn <- cell_of_node[hp]
drr <- abs((cp - 1L) %/% nc - (cn - 1L) %/% nc)
dcc <- abs((cp - 1L) %%  nc - (cn - 1L) %%  nc)
len[hp] <- HIKE_RES_M * sqrt(drr^2 + dcc^2)

# Pointer jumping: follow parents up to the access point, summing hiking distance
acc <- len
repeat {
  nxt <- up[up]
  if (identical(nxt, up)) break
  acc <- acc + acc[up]
  up  <- nxt
}

drive_t_root <- rep(NA_real_, S); drive_t_root[src_nodes] <- rv$drive_time_s
drive_d_root <- rep(NA_real_, S); drive_d_root[src_nodes] <- rv$drive_dist_m
root <- up[seq_len(nn)]

# `dist` is the objective (see OBJ_* in 00_config.R). At the access cell it equals the
# driving part, so everything beyond that is hiking:  hike objective = OBJ_PER_HIKE_MIN * minutes.
objective    <- dist[seq_len(nn)]; objective[!is.finite(objective)] <- NA
drive_time_s <- drive_t_root[root]
drive_dist_m <- drive_d_root[root]
hike_raw     <- (objective - dist[root]) / OBJ_PER_HIKE_MIN * 60
if (any(hike_raw < -1, na.rm = TRUE)) {
  stop(sum(hike_raw < -1, na.rm = TRUE), " cells have a negative hiking time; ",
       "the access-point assignment is wrong. Please report this.")
}
hike_time_s  <- pmax(0, hike_raw)                          # removes rounding noise only
total_time_s <- drive_time_s + hike_time_s
hike_dist_m  <- acc[seq_len(nn)]
unreach <- is.na(objective) | is.na(drive_time_s)
hike_time_s[unreach] <- NA; hike_dist_m[unreach] <- NA; total_time_s[unreach] <- NA

# ---- 5. Back to rasters at HIKE_RES_M -------------------------------------------------
fill <- function(vals) {
  x <- rep(NA_real_, n); x[pass] <- vals
  vec_to_raster(fr_h, x)
}
stk <- c(fill(total_time_s), fill(drive_time_s), fill(hike_time_s),
         fill(drive_dist_m), fill(hike_dist_m), fill(objective))
names(stk) <- c("total_time_s", "drive_time_s", "hike_time_s",
                "drive_dist_m", "hike_dist_m", "objective")
terra::writeRaster(stk, file.path(CACHE_DIR, "cost_hike_res.tif"), overwrite = TRUE)

# ---- 6. Sample at the 200 m cell centroids ------------------------------------------
# Use the hiking cell under the centroid; if it is water/unreachable, use the
# nearest reachable hiking cell inside the same 200 m cell. cost_offset_m says
# how far from the centroid the sampled point is (0 = exactly the centroid cell).
message("Sampling cost at 200 m cell centroids ...")
cent <- terra::xyFromCell(g$template, seq_len(terra::ncell(g$template)))
m    <- CELL_M / HIKE_RES_M
offs <- (seq_len(m) - (m + 1) / 2) * HIKE_RES_M
Tvec <- terra::values(stk[["objective"]], mat = FALSE)    # finite = reachable

best_d    <- rep(Inf, nrow(cent))
best_cell <- rep(NA_integer_, nrow(cent))
any_pass  <- rep(FALSE, nrow(cent))                        # any passable ground in the cell?
for (ox in offs) for (oy in offs) {
  cc <- terra::cellFromXY(stk, cbind(cent[, 1] + ox, cent[, 2] + oy))
  ok <- !is.na(cc)
  any_pass[ok] <- any_pass[ok] | pass[cc[ok]]
  ok[ok] <- is.finite(Tvec[cc[ok]])
  d2  <- ox^2 + oy^2
  upd <- ok & d2 < best_d
  best_d[upd]    <- d2
  best_cell[upd] <- cc[upd]
}
V <- terra::values(stk)
cost <- data.frame(
  cell_id        = seq_len(nrow(cent)),
  drive_time_min = V[best_cell, "drive_time_s"] / 60,
  drive_dist_km  = V[best_cell, "drive_dist_m"] / 1000,
  hike_time_min  = V[best_cell, "hike_time_s"]  / 60,
  hike_dist_km   = V[best_cell, "hike_dist_m"]  / 1000,
  total_time_min = V[best_cell, "total_time_s"] / 60,
  access_objective = V[best_cell, "objective"],
  cost_offset_m  = ifelse(is.finite(best_d), sqrt(best_d), NA_real_)
)
cost$total_dist_km <- cost$drive_dist_km + cost$hike_dist_km
# Why does a cell have no value?
#   no_passable_ground   : every hiking cell inside it is water / impassable / no data
#   unreachable_from_road: has walkable ground, but it is cut off from every road
#                          (e.g. an island, or land enclosed by water / impassable marsh)
cost$cost_status <- ifelse(!is.na(best_cell), "ok",
                           ifelse(any_pass, "unreachable_from_road", "no_passable_ground"))
saveRDS(cost, file.path(CACHE_DIR, "cost_200m.rds"))
print(table(cost$cost_status))

message("Cells with a cost value: ", sum(!is.na(cost$total_time_min)), " of ", nrow(cost))
message("Total time (min): median ", round(median(cost$total_time_min, na.rm = TRUE), 1),
        ", max ", round(max(cost$total_time_min, na.rm = TRUE), 1))
message("Done (", round(difftime(Sys.time(), t0, units = "mins"), 1), " min)")
