# simv2_10_survey_dataset.R
# One realisation of a survey: 100 random land cells, 5 visits each, consecutive visits
# gap_days apart (6 by default). Saves an R object (.rds) and CSV files.
# Run after the world exists:  cfg <- sim_config(...); world <- sim_init(cfg)

# ---- settings ----------------------------------------------------------------
survey_seed <- 123          # seed for the choice of cells and start days
n_cells     <- 100
n_visits    <- 5            # visits per cell
gap_days    <- 6            # days between consecutive visits to the same cell
prefix      <- "sim_survey_100S_5V"

world <- sim_get_world()
cfg   <- world$cfg
sp    <- world$species
ok    <- !world$land$water

out_dir <- dirname(cfg$covariate_file)       # saved next to the grid CSV; change if you like

# ---- grid table: cell_id and travel cost / time -------------------------------
grid <- read.csv(cfg$covariate_file, stringsAsFactors = FALSE)
cost_cols <- c("drive_time_min", "drive_dist_km", "hike_time_min", "hike_dist_km",
               "total_time_min", "total_dist_km", "total_cost_eur")
stopifnot(all(c("cell_id", "row", "col", cost_cols) %in% names(grid)))
grid_key <- paste(grid$row, grid$col)

# ---- 1. random land cells ---------------------------------------------------------
set.seed(survey_seed)
land_cells <- which(ok, arr.ind = TRUE)                       # columns: row, col
pick <- land_cells[sample(nrow(land_cells), n_cells), , drop = FALSE]

cells <- data.frame(site_index = seq_len(n_cells), row = pick[, "row"], col = pick[, "col"])
cells$cell_id <- grid$cell_id[match(paste(cells$row, cells$col), grid_key)]

# standardised covariates of the sampled cells (the values the simulator used)
z_names <- c(forest = "forest_prop_z", edge = "edge_density_z", buffer = "forest_buf1km_z")
for (nm in cfg$cov_names)
  cells[[z_names[[nm]]]] <- world$land$x[cbind(cells$row, cells$col, match(nm, cfg$cov_names))]

# total travel time, distance and cost from the grid pipeline
cells[cost_cols] <- grid[match(cells$cell_id, grid$cell_id), cost_cols]

# ---- 2. visit days: random start, then every gap_days -------------------------------
span  <- gap_days * (n_visits - 1)                            # 24 days for 5 visits and gap_days = 6
start <- sample(cfg$doy_min:(cfg$doy_max - span), n_cells, replace = TRUE)

visits <- data.frame(site_index = rep(seq_len(n_cells), each = n_visits),
                     visit      = rep(seq_len(n_visits), times = n_cells))
visits$cell_id <- cells$cell_id[visits$site_index]
visits$row     <- cells$row[visits$site_index]
visits$col     <- cells$col[visits$site_index]
visits$doy     <- start[visits$site_index] + (visits$visit - 1) * gap_days

stopifnot(max(visits$doy) <= cfg$doy_max,
          all(unlist(tapply(visits$doy, visits$site_index, diff)) == gap_days))

# ---- 3. simulate the recordings (what the sampler would see) --------------------------
det <- sim_collect(data.frame(i = visits$row, j = visits$col, doy = visits$doy))
vid <- rep(seq_len(nrow(visits)), each = nrow(sp))             # sim_collect keeps the visit order
det$cell_id <- visits$cell_id[vid]
det$visit   <- visits$visit[vid]
names(det)[names(det) == "i"] <- "row"
names(det)[names(det) == "j"] <- "col"

# flat table: detections + cell covariates + travel cost
add_cols <- c(unname(z_names), cost_cols)
m <- match(det$cell_id, cells$cell_id)
# detections only: covariates and travel cost go in the metadata sheet (cells)
detections <- det[, c("cell_id", "row", "col", "visit", "doy", "species_id", "max_score", "n_snip_det")]
rownames(detections) <- NULL

# ---- 4. truth at the sampled cells (NOT for the sampler) ---------------------------------
Z <- sapply(seq_len(nrow(sp)), function(k) world$z[, , k][cbind(cells$row, cells$col)])
truth_z <- data.frame(cell_id    = rep(cells$cell_id, times = nrow(sp)),
                      species_id = rep(sp$species_id, each = n_cells),
                      z          = as.vector(Z))

# ---- 5. save ---------------------------------------------------------------------------------
survey <- list(
  detections = detections,   # long: one row per cell x visit x species (observable data)
  cells      = cells,        # one row per sampled cell: ids, standardised covariates, travel cost
  visits     = visits,       # one row per cell x visit: day of year
  truth      = list(z = truth_z, species = sp),
  meta       = list(survey_seed = survey_seed, world_seed = cfg$seed, n_cells = n_cells,
                    n_visits = n_visits, gap_days = gap_days,
                    snippet_interval_min = cfg$snippets$interval_min,
                    snippet_length_min = cfg$snippets$length_min,
                    n_snippets_per_day = sim_n_snippets(cfg), created = Sys.time())
)
saveRDS(survey, file.path(out_dir, paste0(prefix, ".rds")))
write.csv(detections, file.path(out_dir, paste0(prefix, "_detections.csv")), row.names = FALSE)
write.csv(cells, file.path(out_dir, paste0(prefix, "_metadata.csv")), row.names = FALSE)
write.csv(truth_z,    file.path(out_dir, paste0(prefix, "_truth_z.csv")),    row.names = FALSE)
sim_write_species_table(sp, file.path(out_dir, paste0(prefix, "_truth_species.csv")))

cat("Saved to", normalizePath(out_dir), "\n")
cat("detections:", nrow(detections), "rows (", n_cells, "cells x", n_visits, "visits x", nrow(sp), "species )\n")
cat("sampling days range from", min(visits$doy), "to", max(visits$doy), "\n")
head(detections)
head(cells)
# Load later with: survey <- readRDS("<path>/sim_survey_100S_5V.rds")
