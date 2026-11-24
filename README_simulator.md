# Adaptive Sampling: bird detection simulator (R)

Status: implemented in R (base R only). This README was written on 2026-11-20 from the final code
(`simv2_00_config.R`, `simv2_01_create_world.R`, `simv2_02_simulate_validate.R`, `simv2_03_demo.R`) and
replaces `acoustic_simulator_design_v2.md` and `acoustic_simulator_design.md` (v1). The v1 text is kept
in Appendix A. Companion to `grid_pipeline_design.md`.

## 1. Purpose

Simulate acoustic classifier output for a bird community on the 200 m grid, so that adaptive sampling
designs can be tested. The project asks how forest fragmentation affects the species of a community and
uses adaptive sampling to learn it. The simulator therefore keeps the true effects, true occupancy and
true call sources, so any method can be checked against the truth. The sampler only sees what a recorder
and classifier would give it: the highest classifier score per species and the number of 1-minute
snippets with a detection. The dataset made by the survey script is meant to go into an HMSC model; the
model's parameter estimates will be added to the Rmd workflow.

## 2. Files and how to run

| File | Content |
|---|---|
| `simv2_00_config.R` | `sim_config()` (every number of the design), seeding helpers, random-field helpers |
| `simv2_01_create_world.R` | Landscape loading, species table, occupancy, phenology, observation model for one recording day |
| `simv2_02_simulate_validate.R` | Public interface (`sim_init()`, `simulate()`, `sim_collect()`, truth accessors, plots) and `sim_validate()` |
| `simv2_03_demo.R` | Example session |
| Species workflow Rmd | Maps, community, three example species, simulator demo (last delivered as `simv2_09_species_workflow_v3.Rmd`) |
| Survey dataset script | One realised survey saved as `.rds` and CSV files (last delivered as `simv2_10_survey_dataset.R`) |

Run from the folder that holds the files:

```r
for (f in sort(list.files(".", pattern = "^simv2_0[0-2].*\\.R$"))) source(f)   # 00, 01, 02
cfg   <- sim_config(seed = 1, covariate_file = "path/to/grid_200m_all.csv")    # override any value
world <- sim_init(cfg)                    # builds landscape, species, occupancy; kept as the default world
simulate(120, 130, 140)                   # one recording day at row 120, column 130, day 140
v <- sim_validate()                       # checks
```

`sim_config()` merges nested overrides, for example `sim_config(occupancy = list(p_unaffected = 0))`.
Re-sourcing the files keeps an existing world (it lives in `.sim_env`), so call `sim_init(cfg)` again
after changing the config and make a fresh `cfg` object. The Rmd sources the files and builds a new
world every time it is knitted. The survey script needs a world that already exists in the R session.

## 3. Interface

- `sim_init(cfg, species = NULL, set_default = TRUE)` builds the world once: landscape, species table,
  spatial fields, occupancy probabilities `psi` and true states `z`. Pass a species table (see
  `sim_read_species_table()`) to use a hand-edited community.
- `simulate(i, j, doy)` returns one row per species with `species_id`, `max_score` (highest classifier
  score over the day's snippets, 0 = nothing detected) and `n_snip_det` (number of snippets with at least
  one detection of any source, whatever the score). `i` is the row (1 = north), `j` the column
  (1 = west), both 1-based; `doy` runs from 110 to 209 (about 20 April to 28 July). Water cells cannot be
  sampled (error).
- `sim_collect(design)` runs many recordings; `design` is a data frame with columns `i`, `j`, `doy`.
- Truth, never given to the sampler:
  - `sim_truth()`: config, species table, `psi`, `z`, fields, water mask and covariates.
  - `sim_truth_day(i, j, doy)`: for every species, the full detail of one recording: `max_score`,
    `n_snip_det`, `n_true` (true calls in the snippets), `n_true_snips`, `n_vis_snips`, `n_fp_snips`,
    `visitor`, `fp`, `source` of the highest score (`none`, `true`, `visitor`, `fp`), `z`, `act`
    (the species' activity in this cell on this day, 0 to 1).
  - `sim_truth_maps(world, species_id)`: true richness and mean `psi` maps, and for one species its `psi`,
    `z` and field. `sim_plot_truth()` and `sim_image()` draw them (north up, west left).
- Reproducible and stateless: the random numbers of a recording depend only on (seed, i, j, doy), so a
  repeated or reordered query gives an identical answer. The caller's global `.Random.seed` is restored.
- `simulate()` masks `stats::simulate`; call `stats::simulate()` if you need that one.

## 4. Landscape

- 250 x 250 cells of 200 m (50 x 50 km), from the grid pipeline CSV (`covariate_file`), centred on the
  Ylistö campus, Jyväskylä (see `grid_pipeline_design.md`).
- Required columns: `i`, `j` (or `row`, `col`), `forest_prop`, `edge_density` (or `edge_density_m_ha`),
  `forest_buf1km`, `exclude_water`.
- Three covariates, in this order: `forest` (forest cover of the cell), `edge` (forest edge density, forest
  against all non-forest, shorelines included), `buffer` (forest cover within 1 km). All are z-scored again on
  load, to mean 0 and SD 1 over the non-water cells.
- Water: cells with `exclude_water` (water share above 0.20 in the pipeline), cells missing from the file and
  cells with a missing covariate are excluded and never sampled. Their `psi` is 0.

## 5. Community

### 5.1 Species and their traits

50 occupying species. Each trait is drawn once when the world is built (seed stream 2 unless noted):

| Trait | Default | How it is set |
|---|---|---|
| Migration strategy | 17 resident, 17 short-distance, 16 long-distance | Fixed counts, assigned by position, then rows are shuffled, so strategy is independent of everything else |
| Target prevalence | Beta(3, 4.5), clipped to 0.01-0.9 (mean 0.4, SD about 0.18) | One draw per species |
| Profile | rare, intermediate, common | The class of the drawn prevalence: rare up to 0.2, intermediate 0.2-0.5, common above 0.5. With the Beta this gives about 12 %, 60 % and 28 % of the species, so about 6, 30 and 14 of 50, and the counts vary between worlds |
| Confusable | 9 species | Shared out over the profiles in proportion to their size (largest remainder, about 1, 5 and 3), then random species within each profile. Confusable species have a higher false positive rate |
| Affected by forest fragmentation | 80 % affected, 20 % not | Bernoulli draw per species with probability `p_unaffected` of NOT being affected (own seed stream 4). Unaffected species get all three betas equal to 0 together. About 10 of 50 species (95 % range 5 to 16); more than 25 is essentially impossible |
| Covariate effects (betas) | Normal(mean, 0.4) per covariate: forest +0.6, edge -0.4, buffer +0.5 | Effect on logit(occupancy) per 1 SD of the covariate. Drawn independently per covariate, with no truncation, so a few affected species have weak or opposite-sign effects |
| Spatial field scale | none, 300 m, 1 km, 3 km, 10 km | Drawn independently with probability 0.2 each (own seed stream 3); `field_balanced = TRUE` gives near-equal counts instead |
| Phenology | onset, peak, cessation by strategy (section 7) | Uniform draws within the strategy ranges |
| Calling rate `lambda` | log-uniform 30 to 1000 vocalisations per 24 h at full activity | Same range for every profile, so calling rate is independent of prevalence |
| False positive probability `p_fp` | per snippet: 1.4e-4 to 3.6e-4 (confusable: 3.6e-4 to 2e-3) | Uniform within the range |
| Visitor probability `p_visit` | Log-normal, median 0.01, sdlog 0.9, clipped to 0.001-0.05, per day | One draw per species |

Visitor-only species (never occupy a cell) are switched off (`n_visitor_only = 0`). Setting it above 0 adds
that many extra species with profile `visitor_only`, `lambda` 30-1000, higher `p_visit` (median 0.05,
limits 0.02-0.15) and no betas.

### 5.2 The species table

One editable table (`sim_write_species_table()` and `sim_read_species_table()` export and re-import CSV;
rebuild a world with `sim_init(cfg, species = ...)`). Columns:

`species_id` (`sp01` ... `sp50`), `strategy`, `profile`, `confusable`, `frag_affected` (TRUE or FALSE, NA for
visitor-only), `visitor_only`, `prevalence_target`, `beta0` (NA = solved when the world is built),
`beta_forest`, `beta_edge`, `beta_buffer`, `field_scale_m` (NA = no field), `onset`, `peak`, `cessation`,
`shape_kappa`, `lambda`, `p_fp`, `p_visit`.

## 6. Occupancy

- `logit(psi_sc) = beta0_s + sum_k beta_sk x_ck + field_s(c)`, and the latent state is `z_sc ~ Bernoulli(psi_sc)`
  for species `s` in cell `c`. This single state stands for occupancy or breeding; there are no individuals.
- Intercept: `beta0_s` is solved numerically so that the mean `psi` over land cells equals the species'
  target prevalence, so prevalence does not depend on the betas or the field.
- Spatial field: a Gaussian random field with exponential correlation `exp(-d / scale)`, where the scale is the
  e-folding distance of the species (effective range about 3 x scale). With 200 m cells, a 300 m field has a
  neighbour correlation of 0.51 and almost none beyond 1 km, a 10 km field is close to a regional trend.
  Exact simulation by circulant embedding (FFT on a grid padded to twice its size); a warning is given if
  more than 1 % of the spectrum had to be clipped. One FFT set-up is kept per distinct scale.
- Field strength: every field has SD `field_sd = 1` on the logit scale. With `field_standardise = TRUE`
  (default) each field is rescaled to a mean of 0 and an SD of exactly `field_sd` over the land cells, so
  the SD means the same at 300 m and at 10 km. Species with no field have `psi` driven only by the
  covariates.
- Unaffected species have no covariate effect, so their `psi` comes from the intercept and the field only. A
  naive regression can still show an apparent covariate effect for them when their field is large-scale and
  correlated with the covariates; this is intended, as a hard negative control.

## 7. Phenology

- Each species has `onset`, `peak` and `cessation` (day of year) and an activity curve `a(d)`: 0 outside
  `[onset, cessation]`, 1 at the peak, a Beta-shaped bump in between (concentration `kappa = 4`; a larger
  kappa makes it narrower).
- Strategy ranges (uniform draws; peak must be at least 15 days after onset and cessation at least 15 days
  after peak, otherwise the draw is repeated):

| Strategy | Onset | Peak | Cessation |
|---|---|---|---|
| Resident | 85-110 | 125-145 | 185-215 |
| Short-distance | 112-128 | 138-155 | 185-210 |
| Long-distance | 133-150 | 153-170 | 190-215 |

- Cell timing: every species gets a shift per cell, `delta ~ Normal(0, 4 days)`, and uses `a(doy - delta)` in
  that cell, so the same species is early in some cells and late in others. Visitor events use the unshifted
  curve.
- Outside its season a species has `a = 0`: no true calls and no visitor calls, but false positives can still
  happen.

## 8. Observation model for one recording day

### 8.1 Recording schedule

The recorder takes fixed 1-minute snippets, one every 10 minutes, starting at minute 0 of the day: 144
snippets per day (`snippets$interval_min = 15` gives 96). The schedule is the same for every cell and day,
and calling rates are flat over the 24 hours (no diel pattern). The recorder therefore hears 10 % of the
day, and every snippet has the same chance of a call.

### 8.2 Three independent sources per species, drawn snippet by snippet

Let `act` be the species' activity today, `S = 144` and `len = 1` minute.

1. True calls (only where `z = 1`): the number of calls in each snippet is Poisson with mean
   `lambda x act(doy - delta) x len / 1440`. `lambda` is the number of detectable vocalisations (songs or
   calls) in the whole 24 hours; the snippets catch about 10 % of them.
2. Visitor (only where `z = 0`, so never in an occupied cell): the species visits today with probability
   `p_visit x act(doy)`. While it visits it calls at `lambda x act(doy) x visit_frac` calls per day, and these
   are sampled by the snippets in the same way. `visit_frac = 0.021` stands for a visit of about 30 minutes.
3. False positive (always, also outside the season and in cells the species does not occupy): each snippet
   has a false positive with probability `p_fp`.

### 8.3 Scores

- Each true call scores Uniform(0.8, 1); each visitor call Uniform(0.6, 1); each false positive Uniform(0.3, 1).
- The maximum of N such scores on `[lo, hi]` is drawn directly as `lo + (hi - lo) U^(1/N)`.
- `max_score` is the highest score over all sources and snippets (0 = nothing detected). `n_snip_det` counts
  the snippets in which any source produced a detection, whatever its score; it is above 0 exactly when
  `max_score` is above 0.
- About 29 % of false positives score 0.8 or more (0.2 / 0.7), so a high score alone does not rule out a
  false positive.

### 8.4 What the default settings imply

- Daily chance of hearing an occupying species in at least one snippet at full activity:
  `1 - exp(-0.1 x lambda x act)`. With `lambda` from 30 to 1000 this is 95 % at `lambda = 30` and
  essentially 100 % above 100, so every occupied species is heard on nearly every peak-season day. What
  varies between species and cells is the season (activity near onset and cessation is low) and the cell
  timing shift. With a 15-minute schedule the chance at `lambda = 30` is 86 %.
- To target a daily detection probability `p` with the 10-minute schedule, use `lambda = -10 ln(1 - p)`.
- False positives, per day over 144 snippets: 2 to 5 % for normal species and 5 to 25 % for confusable ones
  (`1 - (1 - p_fp)^144`).
- Visitors: the chance of hearing a visiting species is `1 - exp(-0.1 x lambda x act x visit_frac)`. For a
  30-minute visit that is 6 % at `lambda = 30`, 47 % at 300 and 88 % at 1000, 37 % on average over the
  range. With the mean `p_visit` of about 0.014, a species is heard as a visitor on about 0.5 % of the
  cell-days where it is absent.

## 9. Validation (`sim_validate()`)

`sim_validate(world = NULL, n_cells = 200, n_days = 20, plots = TRUE, seed = 1)` prints tables and returns
them in a list:

1. Realised against target prevalence.
2. Covariate effects of the affected species: drawn betas against the config, and recovery by a logistic
   regression of `z` on the covariates (slopes below 1 are expected because the fields add noise).
   2b. The same regression for the unaffected species (true betas 0): mean and maximum absolute estimate
   and the share above 0.2.
3. Species per field scale, and the empirical field autocorrelation over land cell pairs against
   `exp(-d / scale)`, one species per scale, to 10 km (single fields on a finite grid fall below the curve at
   large scales).
4. Mean onset, peak and cessation by strategy.
5. Simulated recordings (random cells and days): (a) detection of occupied cells near the peak by profile,
   with the mean number of snippets with a detection; (b) false positive, visitor and any-detection rates
   where the species is absent; (c) the source of scores of 0.8 or more; (d) calibration of the score (how
   often a score is a real call), and (d2) a check that `n_snip_det > 0` exactly when `max_score > 0`;
   (e) visitor-only species, if any. Groups are the profile, with confusable species marked `+conf`.
6. Reproducibility: the same query gives the same answer after other queries, and the caller's random
   number state is untouched.

Plots: activity curves by strategy, the true richness map and the histogram of maximum scores.

## 10. Companion files and downstream use

- **Rmd workflow.** Header parameters: `sim_dir`, `covariate_file`, `world_seed`, `species_forest`,
  `species_edge`, `species_confusable`, `cost_surface_png`, `cell_seed`, `example_cell`, `one_cell`,
  `one_doy`, `demo_cells`, `demo_doy`. It shows the standardised covariates and the cost surface, the
  community tables, phenology, prevalence, the correlation against distance of the spatial fields, three
  example species (a common forest species, an edge species, a confusable species) and a demonstration of
  `simulate()`, `sim_truth_day()` and `sim_collect()`. By default the demo cells run from species-poor to species-rich, and the occupied example cells are species-rich.
- **Survey dataset script.** One realised survey: 100 random land cells, 5 visits each, 6 days between visits,
  random start days (survey seed 123). It saves `.rds` and CSV files next to the grid CSV: detections
  (`cell_id`, `row`, `col`, `visit`, `doy`, `species_id`, `max_score`, `n_snip_det`), cell metadata
  (standardised covariates and travel time, distance and cost from the grid pipeline), visits, and the truth
  (`z` per cell and species, and the species table, not for the sampler).
- **HMSC (planned).** The survey dataset will be fitted with an HMSC model and its parameter estimates added
  to the Rmd, to compare them with the true betas, the affected or unaffected status and the true `z`. The
  simulator gives scores and snippet counts, not presences, so the rule that turns them into the response
  variable (for example a score threshold) has to be chosen downstream.

## 11. Random-number streams

All streams are derived with `sim_seed()` from the world seed `cfg$seed` (default 2026) and restore the
caller's state:

| Stream | Used for |
|---|---|
| (seed, 2) | species traits (strategy, prevalence, betas, phenology, `lambda`, `p_fp`, `p_visit`, confusable, row order) |
| (seed, 3) | spatial field scale of each species |
| (seed, 4) | affected or unaffected by fragmentation |
| (seed, 11, s) | spatial field of species `s` |
| (seed, 12, s) | occupancy draw `z` of species `s` |
| (seed, 13, i, j) | phenology shifts of all species in cell (i, j) |
| (seed, 20, i, j, doy) | all draws of one recording: visitor events, snippet counts, false positives, scores |

Because the fragmentation and field-scale draws have their own streams, `p_unaffected = 0` (switching the
fragmentation draw off) or a different scale set leaves the other species traits unchanged.

## 12. Configuration (`sim_config()`)

| Parameter | Default | Meaning |
|---|---|---|
| `seed` | 2026 | World seed |
| `nrow`, `ncol`, `cell_m` | 250, 250, 200 | Grid size and cell size in metres |
| `doy_min`, `doy_max` | 110, 209 | Detection days |
| `covariate_file` | `outputs/grid_200m_all.csv` | Grid pipeline CSV |
| `cov_names` | forest, edge, buffer | Covariate names; buffer is forest cover within 1 km |
| `community$n_strategy` | 17, 17, 16 | Resident, short, long-distance species (sum = occupying species) |
| `community$n_visitor_only` | 0 | Extra species that never occupy |
| `community$n_confusable` | 9 | Confusable species in total, shared over profiles in proportion to size |
| `snippets$interval_min`, `length_min` | 10, 1 | One snippet of `length_min` minutes every `interval_min` minutes from minute 0 |
| `occupancy$p_unaffected` | 0.2 | Probability a species is not affected by fragmentation (0 = switch the draw off) |
| `occupancy$beta_mean` | forest 0.6, edge -0.4, buffer 0.5 | Mean effect of the affected species |
| `occupancy$beta_sd` | 0.4 for each | SD of effects across affected species |
| `occupancy$prev_shape`, `prev_limits` | (3, 4.5), (0.01, 0.9) | Beta prevalence and its clipping |
| `occupancy$profile_breaks` | 0.2, 0.5 | Rare up to 0.2, intermediate up to 0.5, common above |
| `occupancy$field_scale_m` | NA, 300, 1000, 3000, 10000 | E-folding distances (NA = no field) |
| `occupancy$field_scale_prob` | 0.2 each | Draw probabilities |
| `occupancy$field_balanced` | FALSE | TRUE = near-equal counts per scale |
| `occupancy$field_sd` | 1 | SD of every field (logit scale) |
| `occupancy$field_standardise` | TRUE | Rescale each field to SD exactly `field_sd` over land |
| `pheno$kappa` | 4 | Narrowness of the activity curve |
| `pheno$cell_shift_sd` | 4 | SD in days of the cell timing shift |
| `pheno$min_gap` | 15 | Minimum days onset to peak and peak to cessation |
| `pheno$ranges` | section 7 | Onset, peak and cessation ranges per strategy |
| `profile$<name>$lambda` | 30 to 1000 (log-uniform) | Calling rate; same for every profile |
| `profile$<name>$p_fp` | 1.4e-4 to 3.6e-4 per snippet | False positive probability |
| `confusable$p_fp` | 3.6e-4 to 2e-3 per snippet | Used instead of the profile range for confusable species |
| `visitor$occupying`, `visitor$visitor_only` | median 0.01 and 0.05, sdlog 0.9 and 0.6, limits 0.001-0.05 and 0.02-0.15 | Daily visitor probability |
| `visitor$visit_frac` | 0.021 | Share of the day a visitor calls at its full rate (about 30 minutes) |
| `obs$true_score`, `visitor_score`, `fp_score` | (0.8, 1), (0.6, 1), (0.3, 1) | Score ranges |

## 13. Changes from the earlier designs

| Topic | v1 | v2 | Now |
|---|---|---|---|
| Community | 20 resident, 25 migrant, 5 visitor; habitat guilds | 45 occupying (15 per strategy) plus 5 visitor-only | 50 occupying (17, 17, 16); visitor-only switched off |
| Species types | Guilds (interior, generalist, edge, open) | Detection profiles abundant, standard, sparse, confusable, with fixed counts | Profiles are prevalence classes (rare, intermediate, common) from the drawn prevalence; confusable is a flag on 9 species |
| Prevalence | Not specified | Beta(2, 6), mean 0.25 | Beta(3, 4.5), mean 0.4 |
| Covariate effects | Guild means, quadratic edge, breeding effects | Normal around community means, every species affected | Same draws, but 20 % of the species get all betas equal to 0 (Bernoulli) |
| Spatial field | 100 m to 1 km range | 30 % none, SD 0.5-1.5, scale 300 m to 5 km | Scale none, 300 m, 1, 3 or 10 km, SD 1, standardised |
| Breeding, individuals, passers-by, flyovers, east-only | Yes | Removed | Not present; one latent state per species and cell |
| Observation | Continuous 24 h, Poisson calls | Continuous 24 h, Poisson calls, Bernoulli false positive and visitor per day | Fixed 1-minute snippets every 10 minutes, flat calling rate, per-snippet false positives, visitors that call for part of a day |
| Calling rate | 1 to 300, median 10 | By profile, 0.4 to 200 | 30 to 1000 vocalisations per 24 h, same for all |
| False positive rate | Poisson, with confusers | 0.01 per day (confusable 0.05-0.15) | Per snippet, 2-5 % per day (confusable 5-25 %) |
| Output | Calls, highest and second-highest score | Highest score | Highest score and number of snippets with a detection |
| Cell timing shift | Not specified | Normal, SD 2.5 days | Normal, SD 4 days |
| Landscape | Dummy | Dummy, later the pipeline | Grid pipeline CSV |
| Code | 7 modules | 9 files | 4 files plus the Rmd and the survey script |

## 14. Known limitations and open items

- Calling is flat over 24 hours. There is no diel pattern, and long days at this latitude made it
  unattractive to build in. Calls arrive independently (Poisson), not in bouts, which overstates how well
  the snippets catch a species for a given number of calls.
- `lambda` is a tuning parameter for detectability, not a literal count of what a bird produces. With 30 to
  1000, occupied species are heard almost every peak-season day, so imperfect detection comes mostly from
  the season and the timing shift.
- Visitors are a placeholder: independent each day and cell, with no link to occupied cells, and a fixed
  30-minute visit.
- False positive levels are defaults, not calibrated to a real classifier. If you have false detections per
  hour of audio, they can replace `p_fp` directly.
- Scores are not calibrated probabilities by construction: only false positives score below 0.8. The
  validation reports how often a score is a real call.
- Affected species draw their betas without truncation, so a few have weak or opposite-sign effects.
- Prevalence profile, migration strategy, calling rate and fragmentation status are drawn independently of
  each other.
- Some comments in the config and demo still carry earlier names or numbers: the `lambda` comment quotes
  detection chances for an earlier range, `covariate_file` refers to `simv2_01_landscape.R`, and the demo
  header says `simv2_08_demo.R`. The numbers in this README follow the code.

---

# Appendix A: the original design (v1, `acoustic_simulator_design.md`)

The text below is the first design, kept unchanged for reference (only its heading levels are lowered by two). It was never implemented as written;
see section 13 for how the design changed.

### Adaptive Sampling: acoustic data simulator design (R)

Status: design agreed in discussion 2026-11-15, no code written yet. Companion to `grid_pipeline_design.md`.

#### Purpose
Simulate acoustic-recognition data for a bird community on the 200 m grid, so that adaptive sampling designs can be tested. The estimand is the effect of forest fragmentation (3 covariates) on occupancy and breeding in the community. The simulator therefore stores the true effects for validation.

#### Interface
- `simulate(i, j, doy)` returns ONE realised draw per species (50 rows) for one 24 h recording at cell (i, j) on day `doy`.
- i = row (1 = north), j = column (1 = west), both 1-based. East = high j.
- Season: 100 detection days, DOY 110 to 209 (about 20 Apr to 28 Jul).
- Randomness is seeded by (seed, i, j, doy, species): stateless and reproducible, so repeated queries return identical values. Nothing is stored as a full cube.
- Latent truth (psi, activity, true states, call classes) is NOT part of the output. It is available through a separate truth accessor so the sampling algorithm cannot use it.

#### Landscape
- 250 x 250 cells (200 m), 50 x 50 km, matching the grid pipeline.
- Water: dummy lake-like blobs, about 10-15 % of cells, excluded from sampling (mirrors the `exclude_water` flag).
- Three dummy covariates for now (later: forest_prop, edge density, forest cover in a 1 km buffer): spatially smooth, moderately correlated (buffer and forest cover positive, edge density hump-shaped against forest cover). All standardised to mean 0, SD 1 over non-water cells.
- A single `load_covariates()` function is the only source of the landscape, so the real pipeline output can be swapped in later.

#### Community: 50 species, one editable species table
- Strategy: 20 resident, 25 migrant, 5 visitor (flyover-only). Counts are editable via a helper that generates the default table.
- Habitat guilds for the 45 non-visitor species: 15 forest interior, 15 forest generalist, 10 edge, 5 open/matrix (landscape is very forested). Visitors have no habitat preference. Strategy and guild are independent.
- Species table columns (all editable): strategy, guild, arrival date and departure date (migrants), 3 covariate effects for occupancy, 3 for breeding, optional quadratic edge term, spatial range, east-only flag, call rate lambda, false-positive rate lambda_fp, confuser species and strength kappa, false-positive score shape.

#### Occupancy
- logit psi = b0 + sum_k b_k x_k (+ quadratic edge for edge species) + Gaussian random field.
- Effects are species-specific, drawn around guild means: interior strongly + for forest and buffer cover and - for edge; generalist mildly +; edge species hump-shaped on edge density (linear + and quadratic -); open species - for forest. "Large-scale forest" preference = the 1 km buffer covariate.
- Placeholder guild means and spread (to be edited): effects of about +/-0.5 to 0.8, SD 0.4 across species.
- Spatial autocorrelation: exponential covariance, species-specific range 100 m to 1 km (note: 100 m is below the 200 m cell size, so effectively independent for those species).
- East-only species: a logistic step on column j multiplies out occupancy in the west.
- Stored truth: species effects, true community-mean effects and true richness response to each covariate.

#### Individuals and phenology
- Every occupied cell holds 2 individuals (switch to set 1). The second is only a device for variability in arrival timing, nothing else depends on the count.
- Individual 1 arrives on the species' arrival date (plus noise); individual 2 arrives 1 to 14 days later (uniform). Residents are present all season.
- A cell is present on a day if at least one individual is present. Call rate is not scaled by the number of individuals.
- Breeding is a Bernoulli draw per individual: logit p_breed = g0 + sum_k g_k x_k, so breeding also responds to fragmentation (more breeding in less fragmented forest). Breeding effects follow the guild pattern, correlated with the occupancy effects but smaller, plus noise. A cell is breeding if either individual breeds. Individuals are independent.
- A breeding individual stays from arrival to departure. A non-breeding individual is a passer-by, present only 1 to 5 days at a random time in the species' window. Non-breeders sing at the same rate as breeders.
- Singing activity over the season rises after arrival, peaks during breeding, then declines after breeding or departure (smooth curve, editable shape).

#### Flyovers
- The 5 visitor species have only flyovers and no resident population.
- Other species also get a small chance of a flyover in cells where they are not occupied, higher near occupied areas (occupancy smoothed over about 1 km).
- At most ONE flyover day per species-cell per season (or none), so a flyover can never produce two correct detections a week apart. A flyover produces a few calls (small Poisson) that are correct detections, scored as true calls.

#### ML observation model (one 24 h recording)
- Plain Poisson, no extra daily variability.
- True calls: N_true ~ Poisson(lambda_species x activity(doy)) when the species is present and active. lambda is heavy-tailed across species (median about 10, range about 1 to 300, some above 100 calls per day). Each true call gets a score ~ Uniform(0.8, 1.0). Recall is built into lambda.
- False positives: N_fp ~ Poisson(lambda_fp,species) on any day, including when the species is absent or out of season. Each gets a score ~ Uniform(0.3, 1.0). Problem species can get a skewed score shape (more mass near 1) and a higher rate. Defaults: median lambda_fp about 0.02 (about 0.6 % of absent species-days with a false positive >= 0.8); problem species 0.5 to 5.
- Confusers: lambda_fp,A += kappa x (true calls of confuser species B that day).
- The ML score is raw in [0, 1]; scores below 0.3 are discarded. Calibration to a classification probability is applied downstream. The simulator's implied P(true | score) can be computed from the stored call classes (Option A: class-conditional scores, calibration emerges from the relative rates of true and false calls).
- Output per species: number of calls >= 0.3, highest score, second-highest score. The true class (correct or false positive) of each retained call is in the truth accessor. Downstream analyses can use only the top score or the top two.
- A false positive can score >= 0.8 with probability about 29 %, so a high score alone does not remove false positives.

#### Breeding confirmation (derived, not simulated)
- Helper `breeding_confirmed()`: a cell-species is confirmed breeding if it has correct (true-call) detections on at least two days separated by at least 7 days.
- Passer-bys and flyovers give single-day or short-stay correct detections and fail this rule by design.
- Whether breeders are confirmed depends on how often the design visits the cell.

#### Suggested module layout (R)
`00_sim_config.R`, `01_landscape.R` (dummy + `load_covariates()`), `02_species_table.R`, `03_latent_fields.R` (psi, breeding, flyover chance), `04_phenology.R`, `05_observe.R` (calls, scores, false positives), `06_simulate.R` (`simulate()`, truth accessor, `breeding_confirmed()`), plus a validation script (prevalence, richness maps, calibration curve, confirmed breeding versus truth).

#### Open items
- Placeholder numbers (guild effect means, SDs, lambda and lambda_fp distributions, activity curve shape, flyover rates) to be reviewed once we see the first simulated maps.
- Replace dummy covariates with forest_prop, edge density and 1 km buffer forest cover from the pipeline.
