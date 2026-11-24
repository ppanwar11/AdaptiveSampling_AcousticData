# simv2_00_config.R
# Bird detection simulator v2: configuration and helpers (base R only).
# Design: project doc "acoustic_simulator_design_v2.md".
#
# Usage:
#   cfg   <- sim_config(seed = 1)            # override any value, nested lists are merged
#   world <- sim_init(cfg)                   # builds landscape, species, occupancy
#   simulate(i, j, doy)                      # one 24 h recording: species_id, max_score

# Keep an existing environment when the files are sourced again, so the world built by
# sim_init() survives re-sourcing.
if (!exists(".sim_env", inherits = FALSE) || !is.environment(.sim_env)) {
  .sim_env <- new.env(parent = emptyenv())
}

# ---------------------------------------------------------------------------
# Configuration. Every number of the design lives here. ####
# Override with e.g. sim_config(occupancy = list(beta_mean = c(forest = 1, edge = 0, buffer = 0)))
# ---------------------------------------------------------------------------
sim_config <- function(...) {
  cfg <- list(
    seed = 2026L,
    nrow = 250L, ncol = 250L, cell_m = 200,   # 250 x 250 cells of 200 m
    doy_min = 110L, doy_max = 209L,           # 100 detection days
    covariate_file = "../outputs/grid_200m_all.csv",  # raw grid pipeline file (row, col, forest_prop,
                                           # edge_density_m_ha, forest_buf1km, exclude_water);
                                           # standardised in simv2_01_landscape.R
    cov_names = c("forest", "edge", "buffer"), #buffer is forest cover in 1km radius

    community = list(
      n_strategy = c(resident = 17, short = 17, long = 16),   # occupying species (sum = 50)
      n_visitor_only = 0,                                      # species that never occupy (0 = none; set e.g. 5 to bring them back)
      # profile = prevalence class of the drawn target prevalence (see occupancy$profile_breaks),
      # so the number of species per profile is random, not fixed.
      # confusable species (higher false positive rate): this many in total, shared out over the
      # profiles in proportion to their size (largest remainder), random species within each
      n_confusable = 9
    ),

    # Fixed recording schedule: one snippet of length_min minutes every interval_min minutes,
    # starting at minute 0 of the day (10 min -> 144 snippets, 15 min -> 96 snippets).
    # Calling rates are flat over the 24 h. No diel pattern yet.
    snippets = list(interval_min = 10, length_min = 1),

    occupancy = list(
      # Forest fragmentation: each occupying species is NOT affected with probability p_unaffected
      # (Bernoulli draw); unaffected species get all three betas = 0. 0.2 -> about 10 of 50 species.
      # Set p_unaffected = 0 to switch the draw off (every species gets betas).
      p_unaffected = 0.2,
      # effects of the affected species (per 1 SD of the covariate on the logit scale)
      beta_mean = c(forest = 0.6, edge = -0.4, buffer = 0.5), # community mean effect per covariate
      beta_sd   = c(forest = 0.4, edge = 0.4, buffer = 0.4),  # SD across species per covariate
      # target prevalence of each species ~ Beta(a, b), clipped (mean 0.4, SD about 0.18)
      prev_shape = c(3, 4.5),
      prev_limits = c(0.01, 0.9),
      # profile = prevalence class: rare <= 0.2 < intermediate <= 0.5 < common
      # (with Beta(3, 4.5): about 12% rare, 60% intermediate, 28% common)
      profile_breaks = c(0.2, 0.5),
      # spatial fields. Each occupying species gets one e-folding distance from this set,
      # corr = exp(-d / scale); NA = no spatial field. Practical range (corr ~ 0.05) = 3 x scale.
      field_scale_m = c(NA, 300, 1000, 3000, 10000),
      field_scale_prob = rep(0.2, 5),                         # draw probabilities (independent per species)
      field_balanced = FALSE,                                 # TRUE = near-equal counts per scale, shuffled
      field_sd = 1,                                           # SD of every field on the logit scale
      field_standardise = TRUE                                # rescale each field to SD exactly field_sd over land cells
    ),

    pheno = list(
      kappa = 4,            # width of the Beta-shaped activity curve (larger = narrower)
      cell_shift_sd = 4,  # days, per species and cell. This introduces variation among individuals of same spp
      min_gap = 15,         # minimum days onset -> peak and peak -> cessation
      ranges = list(
        resident = list(onset = c(85, 110),  peak = c(125, 145), cessation = c(185, 215)),
        short    = list(onset = c(112, 128), peak = c(138, 155), cessation = c(185, 210)),
        long     = list(onset = c(133, 150), peak = c(153, 170), cessation = c(190, 215))
      )
    ),

    # Profiles are prevalence classes (share of land cells occupied), see occupancy$profile_breaks.
    # lambda = expected number of detectable vocalisations (songs or calls) per 24 h at activity 1
    #          (log-uniform range), spread evenly over the day. With one 1-min snippet every
    #          10 min, a species is heard on a given day with probability about 1 - exp(-lambda / 10):
    #          lambda 10 -> 63%, 30 -> 95%, 100 or more -> ~100%.
    #          Use lambda = -10 * log(1 - p) to target a daily detection probability p.
    #          The same range for every profile, so calling rate is independent of prevalence;
    #          give a profile its own range to tie them (e.g. rare = c(10, 100)).
    # p_fp   = false positive probability PER SNIPPET (a range = uniform). 1.4e-4 to 3.6e-4 per
    #          snippet is about 0.02 to 0.05 per day with 144 snippets (1 - (1 - p)^144).
    # Visitor-only species (if any) need a lambda too: it is their calling rate while visiting.
    profile = list(
      rare         = list(lambda = c(30, 1000), p_fp = c(1.4e-4, 3.6e-4)),
      intermediate = list(lambda = c(30, 1000), p_fp = c(1.4e-4, 3.6e-4)),
      common       = list(lambda = c(30, 1000), p_fp = c(1.4e-4, 3.6e-4)),
      visitor_only = list(lambda = c(30, 1000), p_fp = c(1.4e-4, 3.6e-4))
    ),

    # For confusable species (confused with some background sounds other species/ or geophony)
    # use this p_fp range (per snippet) instead of their profile's;
    # 3.6e-4 to 2e-3 per snippet is about 0.05 to 0.25 per day with 144 snippets #0.1/day = 7.3e-4/snippet. 0.3/day = 2.5e-3/snippet
    confusable = list(p_fp = c(3.6e-4, 2e-3)),

    # daily visitor probability (before multiplying by the activity curve): log-normal, clipped.
    # A visitor calls at lambda x activity x visit_frac calls per day (placeholder, to be revisited).
    visitor = list(
      occupying    = list(median = 0.01, sdlog = 0.9, limits = c(0.001, 0.05)),
      visitor_only = list(median = 0.05, sdlog = 0.6, limits = c(0.02, 0.15)),
      visit_frac   = 0.021 ## Assuming a visiting individual in unoccuppied cell can be detected in 30 of 1440snippets
    ),

    # score ranges (classifier probability)
    obs = list(true_score = c(0.8, 1.0), visitor_score = c(0.6, 1.0), fp_score = c(0.3, 1.0))
  )
  cfg <- utils::modifyList(cfg, list(...))
  cfg
}

# ---------------------------------------------------------------------------
# Reproducible random numbers ####
# ---------------------------------------------------------------------------
# Combine integers into one valid RNG seed (no integer overflow, order matters).
sim_seed <- function(...) {
  v <- as.numeric(c(...))
  h <- 12345
  for (x in v) h <- (h * 1000003 + x + 1) %% 2147483629
  as.integer(h) + 1L
}

# Evaluate expr with set.seed(seed), then restore the caller's RNG state.
sim_with_seed <- function(seed, expr) {
  has_old <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_old) old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (has_old) {
      assign(".Random.seed", old, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

sim_loguniform <- function(n, r) exp(runif(n, log(r[1]), log(r[2])))

# ---------------------------------------------------------------------------
# Gaussian random field, unit variance, exponential correlation exp(-d / scale_m). ####
# Exact simulation by circulant embedding (FFT) on a grid padded to 2x in each direction.
# ---------------------------------------------------------------------------
sim_grf_setup <- function(nr, nc, cell_m, scale_m) {
  M <- 2L * nr
  N <- 2L * nc
  a <- 0:(M - 1L)
  b <- 0:(N - 1L)
  dx <- pmin(a, M - a) * cell_m
  dy <- pmin(b, N - b) * cell_m
  D <- sqrt(outer(dx^2, dy^2, "+"))
  C <- exp(-D / scale_m)
  lam <- Re(fft(C))
  neg_frac <- -sum(lam[lam < 0]) / (M * N)   # share of variance lost by clipping (should be ~0)
  lam[lam < 0] <- 0
  list(nr = nr, nc = nc, M = M, N = N, sqrt_lam = sqrt(lam), neg_frac = neg_frac)
}

sim_grf_draw <- function(setup) {
  n <- setup$M * setup$N
  xi <- matrix(complex(real = rnorm(n), imaginary = rnorm(n)), setup$M, setup$N)
  Y <- fft(setup$sqrt_lam * xi) / sqrt(n)
  Re(Y)[seq_len(setup$nr), seq_len(setup$nc)]
}

# Standardise a matrix to mean 0, SD 1 over cells where ok is TRUE; other cells set to 0.
sim_standardise <- function(m, ok) {
  v <- m[ok]
  m <- (m - mean(v)) / sd(v)
  m[!ok] <- 0
  m
}
