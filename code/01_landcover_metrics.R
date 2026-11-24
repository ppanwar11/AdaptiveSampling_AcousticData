# =============================================================================
# 01_landcover_metrics.R
# Forest cover, water, marsh and forest edge density per 200 m cell, plus forest cover
# in a circular buffer (BUF_RADIUS_M) around each cell centre.
# Output: cache/metrics_200m.tif  (250 x 250 cells, 9 layers; see METRIC_NAMES)
# =============================================================================
source("code/00_config.R")
t0 <- Sys.time()

g <- make_grid()
message("Centre (EPSG:3067): ", paste(round(g$centre), collapse = ", "),
        "  -> snapped: ", paste(g$snapped_centre, collapse = ", "))
message("Grid extent: ", paste(as.vector(g$ext), collapse = ", "))

stopifnot(CELL_M %% CLC_RES_M == 0)
fact  <- CELL_M / CLC_RES_M          # 10 Corine pixels per cell side
npix  <- fact^2                       # 100 pixels per cell

# ---- 1. Read Corine for the grid + 1 km buffer ----------------------------------
message("Reading Corine window ...")
win <- read_clc_window(buffer_ext(g$ext, EDGE_BUFFER_M))
v   <- win$codes
nr  <- nrow(win$geom); nc <- ncol(win$geom)

ok        <- !is.na(v)
is_forest <- ok & (v %in% FOREST_CODES)
is_water  <- ok & (v %in% WATER_CODES)
is_marsh  <- ok & (v %in% MARSH_CODES)

# Binary forest map: 1 forest, 0 anything else, NA = no data
Fv  <- ifelse(ok, as.numeric(is_forest), NA_real_)
# Same, but water treated as "no data" so forest/water boundaries are ignored
Fnw <- Fv; Fnw[is_water] <- NA_real_

# ---- 2. Edge length per pixel -----------------------------------------------------
# Edge = a boundary between a forest pixel and a non-forest pixel (4-neighbour
# adjacency).  Every differing pair has length CLC_RES_M; half of it is given to
# each of the two pixels, so edges that run along a cell boundary are shared
# fairly between the two cells.  Pairs involving NA are ignored.
edge_length_per_pixel <- function(x, nr, nc, px) {
  M  <- matrix(x, nrow = nr, ncol = nc, byrow = TRUE)
  e  <- matrix(0, nr, nc)
  dh <- M[, -nc, drop = FALSE] != M[, -1, drop = FALSE]   # horizontal pairs
  dv <- M[-nr, , drop = FALSE] != M[-1, , drop = FALSE]   # vertical pairs
  dh[is.na(dh)] <- FALSE
  dv[is.na(dv)] <- FALSE
  e[, -nc] <- e[, -nc] + dh
  e[, -1]  <- e[, -1]  + dh
  e[-nr, ] <- e[-nr, ] + dv
  e[-1, ]  <- e[-1, ]  + dv
  as.vector(t(e * (px / 2)))                              # back to row-major
}

message("Computing edges ...")
edge_all <- edge_length_per_pixel(Fv,  nr, nc, CLC_RES_M)
edge_nw  <- edge_length_per_pixel(Fnw, nr, nc, CLC_RES_M)

# ---- 3. Aggregate 20 m -> 200 m (crop away the buffer first) -------------------
agg_sum <- function(x) {
  r <- vec_to_raster(win$geom, x)
  r <- terra::crop(r, g$ext)
  terra::aggregate(r, fact = fact, fun = "sum", na.rm = TRUE)
}

forest_n <- agg_sum(as.numeric(is_forest))
water_n  <- agg_sum(as.numeric(is_water))
marsh_n  <- agg_sum(as.numeric(is_marsh))
valid_n  <- agg_sum(as.numeric(ok))
edge_len <- agg_sum(edge_all)
edge_len_nw <- agg_sum(edge_nw)

cell_ha  <- CELL_M^2 / 10000
land_n   <- valid_n - water_n

forest_prop      <- forest_n / npix
forest_prop_land <- terra::ifel(land_n > 0, forest_n / land_n, NA)
edge_density     <- edge_len    / cell_ha
edge_density_nw  <- edge_len_nw / cell_ha
water_prop       <- water_n / npix
marsh_prop       <- marsh_n / npix
valid_prop       <- valid_n / npix

# ---- 3b. Forest cover in a circle around each cell centre -----------------------------
# Counted directly on the 20 m pixels read above (grid + EDGE_BUFFER_M). A pixel is in
# the circle if its centre is within BUF_RADIUS_M of the cell centre. Cell centres lie
# on pixel corners (CELL_M / CLC_RES_M is even), so every cell has the same circle
# pattern: the k-th pixel row above and below the cell's centre line spans h_k pixels
# on each side of its centre column. Row sums come from cumulative sums along each pixel
# row, so there is no loop over cells or pixels (50 steps for a 1 km circle).
message("Computing ", BUF_RADIUS_M, " m buffer ...")
R   <- BUF_RADIUS_M / CLC_RES_M          # radius in pixels
off <- EDGE_BUFFER_M / CLC_RES_M         # pixels between window edge and grid edge
stopifnot(abs(R - round(R)) < 1e-9, abs(off - round(off)) < 1e-9,
          fact %% 2 == 0,
          EDGE_BUFFER_M + CELL_M / 2 >= BUF_RADIUS_M,   # circle must stay inside the window
          nr == N_CELLS * fact + 2 * off, nc == N_CELLS * fact + 2 * off)
R   <- round(R); off <- round(off)
cy  <- as.integer(round(off + fact * (seq_len(N_CELLS) - 0.5)))  # centre line below pixel row cy
cx  <- cy                                                         # same for columns (square grid)
# Half-width (pixels per side) of the circle in the k-th row away from the centre line:
# pixel centres are at (m - 1/2, k - 1/2) pixels, inside if (2m-1)^2 + (2k-1)^2 <= 4 R^2
hk  <- floor((sqrt(4 * R^2 - (2 * seq_len(R) - 1)^2) + 1) / 2)
circle_n <- 4 * sum(hk)                  # pixels in the circle (2 rows x 2 sides per k)

circle_count <- function(x) {            # x: nr x nc matrix of 0/1; returns N_CELLS x N_CELLS
  CS  <- cbind(0, t(apply(x, 1, cumsum)))              # CS[r, c + 1] = sum of x[r, 1..c]
  tot <- matrix(0, N_CELLS, N_CELLS)
  for (k in seq_len(R)) {
    a  <- cx - hk[k] + 1L;  b  <- cx + hk[k]           # columns a..b (by cell column)
    r1 <- cy - k + 1L;      r2 <- cy + k               # pixel rows above / below the line
    tot <- tot + CS[r1, b + 1L] - CS[r1, a] + CS[r2, b + 1L] - CS[r2, a]
  }
  tot
}
to_cells <- function(m) vec_to_raster(forest_n, as.vector(t(m)))     # matrix -> same geometry as forest_n

buf_forest_n <- circle_count(matrix(as.numeric(is_forest), nr, nc, byrow = TRUE))
buf_valid_n  <- circle_count(matrix(as.numeric(ok),        nr, nc, byrow = TRUE))
forest_buf     <- to_cells(buf_forest_n / circle_n)    # water and no-data count as non-forest
buf_valid_prop <- to_cells(buf_valid_n  / circle_n)    # share of the circle with Corine data
message("  circle = ", circle_n, " pixels (", round(circle_n * CLC_RES_M^2 / 1e6, 3),
        " km2; a true circle of this radius is ",
        round(pi * BUF_RADIUS_M^2 / 1e6, 3), " km2)")

metrics <- c(forest_prop, forest_prop_land, edge_density, edge_density_nw,
             water_prop, marsh_prop, valid_prop, forest_buf, buf_valid_prop)
names(metrics) <- METRIC_NAMES

# The aggregated raster must coincide with the study-grid template
stopifnot(terra::compareGeom(metrics, g$template, crs = FALSE,
                             stopOnError = FALSE))

out <- file.path(CACHE_DIR, "metrics_200m.tif")
terra::writeRaster(metrics, out, overwrite = TRUE)

# ---- 4. Quick summary --------------------------------------------------------------
s <- terra::values(metrics)
message("\nSummary of the ", nrow(s), " cells:")
print(round(apply(s, 2, function(z) c(min = min(z, na.rm = TRUE),
                                      mean = mean(z, na.rm = TRUE),
                                      max = max(z, na.rm = TRUE))), 3))
message("Cells with water_prop > ", WATER_EXCLUDE_THRESHOLD, ": ",
        sum(s[, "water_prop"] > WATER_EXCLUDE_THRESHOLD, na.rm = TRUE),
        " of ", nrow(s))
message("Wrote ", out, "  (", round(difftime(Sys.time(), t0, units = "mins"), 1), " min)")

