# Adaptive sampling grid: forest, edge, water and travel cost (Ylistö, Jyväskylä)

50 x 50 km grid of 250 x 250 cells (200 m), centred on Ylistö campus (62.2293 N, 25.7414 E),
EPSG:3067. Land cover: Corine Land Cover 2024, 20 m.

## Run
1. Install: `install.packages(c("terra","sf","igraph","dodgr","osmextract"))`
2. Edit `CLC_PATH` (and anything else you want) in `00_config.R`.
3. `setwd("path/to/pipeline"); source("run_all.R")`
   Steps can also be run one at a time; each saves its result in `cache/`.

| Step | Script | Result in `cache/` |
|---|---|---|
| 1 | `01_landcover_metrics.R` | `metrics_200m.tif` |
| 2 | `02_driving.R` | `drive_vertices.rds` (+ cached OSM roads) |
| 3 | `03_hiking_cost.R` | `cost_200m.rds`, `cost_hike_res.tif` |
| 4 | `04_export.R` | final files in `outputs/` |

## Outputs (`outputs/`)
- `grid_200m_all.csv`: all 62,500 cells, with an `exclude_water` flag
- `grid_200m_filtered.csv`: only cells with water <= 20 %
- `grid_200m.tif`: multi-band GeoTIFF; cells with water > 20 % are NA (`grid_200m_unmasked.tif` too if `WRITE_UNMASKED_TIF`)
- `grid_landcover_overview.png`: forest, edge and water maps + histograms (excluded cells grey)
- `grid_cost_surface.png`: driving time/distance, hiking time/distance, total time, fuel cost

## Columns
| Column | Meaning |
|---|---|
| `cell_id`, `row`, `col` | cell number (row-major from top-left), row, column |
| `x_3067`, `y_3067`, `lon`, `lat` | cell centroid in EPSG:3067 and WGS84 |
| `forest_prop` | share of the whole cell that is forest (0-1) |
| `forest_prop_land` | share of the non-water area that is forest |
| `edge_density_m_ha` | forest / non-forest boundary length per hectare (water shorelines count as edges) |
| `edge_density_nowater_m_ha` | same, ignoring forest-water boundaries |
| `water_prop`, `marsh_prop`, `valid_prop` | open water, marsh/peatbog, and valid-data shares |
| `exclude_water` | TRUE if `water_prop` > 0.20 |
| `drive_time_min`, `drive_dist_km` | car part: campus to the road access point used |
| `hike_time_min`, `hike_dist_km` | walking part: access point to the cell |
| `total_time_min`, `total_dist_km` | drive + hike, minimised over all possible access points |
| `access_objective` | value of the quantity minimised when choosing the access point (see `OBJ_*` in `00_config.R`) |
| `fuel_l`, `fuel_cost_eur` | fuel use and cost from `drive_dist_km`, `FUEL_L_PER_100KM`, `FUEL_PRICE_EUR_PER_L`; doubled if `ROUND_TRIP` |
| `time_cost_eur`, `total_cost_eur` | time cost (`VALUE_OF_TIME_EUR_H`, 0 by default) and fuel + time cost |
| `cost_status` | `ok`, `no_passable_ground` (cell is all water / impassable) or `unreachable_from_road` (walkable but cut off from every road) |
| `cost_offset_m` | distance from the cell centroid to the sampled hiking cell (0 = centroid itself) |

## Class definitions (Corine 2024, LEVEL4)
- Forest: 3111, 3112, 3121, 3122, 3123, 3131, 3132, 3133 (all forest classes incl. peatland and rocky).
  Not forest: transitional woodland/shrub 3241-3246, agro-forestry 2441, urban green 1411.
- Open water: 5111, 5121, 5231.  Marsh/bog: 4111, 4112, 4121, 4211, 4212.
- The legend carries canopy-cover thresholds but no tree height, so the 5 m rule cannot be checked directly.

## Method notes
- Grid is snapped to the 200 m lattice, so each cell holds exactly 10 x 10 Corine pixels.
- Edge density: 4-neighbour pixel pairs with different forest status, 20 m each; edges between two cells are shared half-half. Land cover is read with a 1 km buffer so border cells are correct.
- Driving: OSM roads (car profile, `dodgr`), grid + 10 km buffer. Private roads dropped (`ALLOW_PRIVATE_ROADS`); access "no" always dropped.
  Forest roads (highway=track and forestry/agricultural access) are OFF by default: switch `ALLOW_FOREST_ROADS` in `00_config.R`, then re-run steps 02-04.
- Hiking: 4 km/h x friction on a 40 m graph (8 neighbours); open water impassable; marsh/bog/wet ground walked at `MARSH_SPEED_KMH` (1 km/h; `MARSH_IMPASSABLE` to block it); peatland forest 1.5x.
  Roads that cross open water (bridges, causeways) are walkable and can be access points (`BRIDGES_PASSABLE`).
  Rivers narrower than about one hiking cell may disappear when the 20 m map is aggregated to 40 m.
- Cells whose centroid is water or unreachable take the nearest reachable hiking cell inside the same 200 m cell.

## Known limits
- Corine 2024 is a mixed-source 20 m product: narrow roads, power lines and small clearings are under-represented, so absolute edge densities are lower than from finer data. Rankings across cells are the intended use.
- OSM speed limits are incomplete; dodgr fills gaps with defaults by road class. Digiroad would be better for speeds.
- Travel times assume summer conditions and no terrain slope.
