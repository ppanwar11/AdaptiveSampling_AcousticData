# ===========================================================================
# Landscape: water mask and three standardised covariates from the grid pipeline CSV####
# ===========================================================================
# A landscape is a list:
#   nrow, ncol, cell_m
#   water : logical matrix [nrow, ncol], TRUE = excluded (not sampled)
#   x     : numeric array [nrow, ncol, k], standardised (mean 0, SD 1 over non-water cells)
#   X     : numeric matrix [nrow * ncol, k], same values, one row per cell (column-major cell order)
# i = row (1 = north), j = column (1 = west).
#
# Required CSV columns: i, j, forest_prop, edge_density, forest_buf1km, exclude_water
# The covariates are used in the order of cfg$cov_names (forest, edge, buffer) and are
# z-scored again on load.

sim_load_covariates <- function(cfg) {
  if (is.null(cfg$covariate_file)) stop("set covariate_file in sim_config()")
  d <- utils::read.csv(cfg$covariate_file, stringsAsFactors = FALSE)
  if (!"i" %in% names(d) && "row" %in% names(d)) names(d)[names(d) == "row"] <- "i"
  if (!"j" %in% names(d) && "col" %in% names(d)) names(d)[names(d) == "col"] <- "j"
  if (!"edge_density" %in% names(d) && "edge_density_m_ha" %in% names(d))
    names(d)[names(d) == "edge_density_m_ha"] <- "edge_density"
  need <- c("i", "j", "forest_prop", "edge_density", "forest_buf1km", "exclude_water")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("covariate_file lacks columns: ", paste(miss, collapse = ", "))

  nr <- cfg$nrow
  nc <- cfg$ncol
  if (any(d$i < 1 | d$i > nr | d$j < 1 | d$j > nc)) stop("covariate_file has i or j outside the grid")
  idx <- cbind(d$i, d$j)
  ew <- d$exclude_water
  ew <- if (is.character(ew)) toupper(ew) %in% c("TRUE", "T", "1") else as.logical(ew)
  ew[is.na(ew)] <- TRUE
  water <- matrix(TRUE, nr, nc)          # cells missing from the file are treated as excluded
  water[idx] <- ew
  put <- function(v) { m <- matrix(NA_real_, nr, nc); m[idx] <- v; m }
  covs <- list(put(d$forest_prop), put(d$edge_density), put(d$forest_buf1km))
  water <- water | Reduce(`|`, lapply(covs, is.na))
  ok <- !water
  covs <- lapply(covs, function(m) sim_standardise(m, ok))

  k <- length(covs)
  x <- array(unlist(covs), dim = c(nr, nc, k), dimnames = list(NULL, NULL, cfg$cov_names))
  X <- matrix(x, nr * nc, k, dimnames = list(NULL, cfg$cov_names))
  list(nrow = nr, ncol = nc, cell_m = cfg$cell_m, water = water, x = x, X = X)
}

# ===========================================================================
# The editable species table ####
# ===========================================================================
# Columns:
#   species_id, strategy (resident / short / long),
#   profile (prevalence class: rare / intermediate / common, or visitor_only),
#   confusable (TRUE = higher false positive rate; spread over the profiles in proportion to size), visitor_only,
#   frag_affected (TRUE = affected by forest fragmentation; FALSE = all betas are 0; NA = visitor-only),
#   prevalence_target, beta0 (NA = solved in sim_init), beta_<cov> (one per covariate),
#   field_scale_m (e-folding distance in m, NA = no spatial field; the field SD is cfg$occupancy$field_sd),
#   onset, peak, cessation, shape_kappa,
#   lambda (expected true calls per 24 h at activity 1), p_fp (false positive probability per snippet),
#   p_visit (daily visitor probability in an unoccupied cell for the species).

sim_make_species <- function(cfg, land) {
  sim_with_seed(sim_seed(cfg$seed, 2L), {
    com <- cfg$community
    occ <- cfg$occupancy
    ph <- cfg$pheno
    strat_names <- names(com$n_strategy)
    n_occ <- sum(com$n_strategy)
    n_vo <- com$n_visitor_only
    n <- n_occ + n_vo
    k <- length(cfg$cov_names)
    classes <- c("rare", "intermediate", "common")

    strategy <- c(rep(strat_names, com$n_strategy),
                  sample(strat_names, n_vo, replace = TRUE))
    visitor_only <- c(rep(FALSE, n_occ), rep(TRUE, n_vo))

    # target prevalence ~ Beta, clipped; the profile is its class (rare / intermediate / common)
    prev <- pmin(pmax(rbeta(n_occ, occ$prev_shape[1], occ$prev_shape[2]),
                      occ$prev_limits[1]), occ$prev_limits[2])
    profile <- c(as.character(cut(prev, c(-Inf, occ$profile_breaks, Inf), labels = classes)),
                 rep("visitor_only", n_vo))
    prev <- c(prev, rep(0, n_vo))          # visitor-only species never occupy

    # confusable species: n_confusable in total, shared out over the profiles in proportion to
    # their size (largest remainder), then a random subset within each profile
    confusable <- rep(FALSE, n)
    sizes <- vapply(classes, function(g) as.numeric(sum(profile == g)), numeric(1))
    n_conf <- sim_allocate(sum(com$n_confusable), sizes)
    for (g in classes) {
      if (n_conf[[g]] == 0) next
      idx <- which(profile == g)
      confusable[idx[sample.int(length(idx), n_conf[[g]])]] <- TRUE
    }

    # covariate effects: beta_sk ~ Normal(mu_k, sigma_k); none for visitor-only species
    B <- matrix(NA_real_, n, k)
    for (q in seq_len(k)) {
      nm <- cfg$cov_names[q]
      B[, q] <- rnorm(n, occ$beta_mean[[nm]], occ$beta_sd[[nm]])
    }
    B[visitor_only, ] <- NA_real_

    # forest fragmentation: Bernoulli draw per occupying species, with probability
    # occ$p_unaffected of NOT being affected; unaffected species get all betas = 0.
    # The betas above are drawn for every species first, and the Bernoulli draw has its own seed
    # stream, so p_unaffected = 0 (draw switched off) leaves all other traits unchanged.
    unaffected <- rep(FALSE, n)
    if (n_occ > 0 && isTRUE(occ$p_unaffected > 0))
      unaffected[seq_len(n_occ)] <- sim_with_seed(sim_seed(cfg$seed, 4L),
                                                  runif(n_occ) < occ$p_unaffected)
    B[unaffected, ] <- 0
    frag_affected <- ifelse(visitor_only, NA, !unaffected)

    # e-folding distance of the spatial field: one value per species from occ$field_scale_m
    # (NA = no field). Own seed stream, so changing the scale set does not shift other traits.
    field_scale <- rep(NA_real_, n)
    if (n_occ > 0) {
      sc <- occ$field_scale_m
      pr <- occ$field_scale_prob
      if (length(pr) != length(sc)) stop("field_scale_prob must have one value per field_scale_m")
      pick <- sim_with_seed(sim_seed(cfg$seed, 3L), {
        if (isTRUE(occ$field_balanced)) sample(rep_len(seq_along(sc), n_occ))
        else sample(seq_along(sc), n_occ, replace = TRUE, prob = pr)
      })
      field_scale[seq_len(n_occ)] <- sc[pick]
    }

    P <- matrix(NA_real_, n, 3, dimnames = list(NULL, c("onset", "peak", "cessation")))
    for (s in seq_len(n)) P[s, ] <- sim_draw_pheno(strategy[s], ph)

    lambda <- rep(NA_real_, n)
    p_fp <- rep(NA_real_, n)
    for (g in unique(profile)) {
      idx <- which(profile == g)
      pr <- cfg$profile[[g]]
      if (!anyNA(pr$lambda)) lambda[idx] <- sim_loguniform(length(idx), pr$lambda)
      for (cf in c(FALSE, TRUE)) {       # confusable species use cfg$confusable$p_fp
        j <- idx[confusable[idx] == cf]
        if (!length(j)) next
        rng <- if (cf) cfg$confusable$p_fp else pr$p_fp
        p_fp[j] <- if (length(rng) == 2) runif(length(j), rng[1], rng[2]) else rng
      }
    }

    p_visit <- numeric(n)
    for (vo in c(FALSE, TRUE)) {
      idx <- which(visitor_only == vo)
      pv <- if (vo) cfg$visitor$visitor_only else cfg$visitor$occupying
      p_visit[idx] <- pmin(pmax(rlnorm(length(idx), log(pv$median), pv$sdlog), pv$limits[1]), pv$limits[2])
    }

    sp <- data.frame(species_id = NA_character_, strategy = strategy, profile = profile,
                     confusable = confusable, frag_affected = frag_affected,
                     visitor_only = visitor_only, prevalence_target = prev,
                     beta0 = NA_real_, stringsAsFactors = FALSE)
    for (q in seq_len(k)) sp[[paste0("beta_", cfg$cov_names[q])]] <- B[, q]
    sp$field_scale_m <- field_scale
    sp$onset <- P[, "onset"]
    sp$peak <- P[, "peak"]
    sp$cessation <- P[, "cessation"]
    sp$shape_kappa <- ph$kappa
    sp$lambda <- lambda
    sp$p_fp <- p_fp
    sp$p_visit <- p_visit

    # shuffle the rows so visitor-only species are not recognisable by their id
    sp <- sp[sample.int(n), , drop = FALSE]
    rownames(sp) <- NULL
    sp$species_id <- sprintf("sp%02d", seq_len(n))
    sp
  })
}

# Split an integer total over groups in proportion to their sizes (largest remainder method).
# Named result, never more than the group size.
sim_allocate <- function(total, sizes) {
  total <- min(total, sum(sizes))
  if (total <= 0) return(sizes * 0)
  q <- total * sizes / sum(sizes)
  out <- floor(q)
  rest <- total - sum(out)
  if (rest > 0) {
    add <- order(q - out, decreasing = TRUE)[seq_len(rest)]
    out[add] <- out[add] + 1
  }
  out
}

# onset < peak < cessation inside the strategy ranges, whole days, minimum gaps.
sim_draw_pheno <- function(strategy, ph) {
  r <- ph$ranges[[strategy]]
  for (try in 1:2000) {
    onset <- round(runif(1, r$onset[1], r$onset[2]))
    peak <- round(runif(1, r$peak[1], r$peak[2]))
    cess <- round(runif(1, r$cessation[1], r$cessation[2]))
    if (peak - onset >= ph$min_gap && cess - peak >= ph$min_gap)
      return(c(onset = onset, peak = peak, cessation = cess))
  }
  stop("could not draw onset < peak < cessation for strategy ", strategy,
       "; check pheno$ranges and min_gap")
}

sim_write_species_table <- function(sp, file) {
  utils::write.csv(sp, file, row.names = FALSE)
  invisible(file)
}

sim_read_species_table <- function(file) {
  sp <- utils::read.csv(file, stringsAsFactors = FALSE)
  need <- c("species_id", "strategy", "profile", "visitor_only", "prevalence_target",
            "field_scale_m", "onset", "peak", "cessation", "shape_kappa",
            "lambda", "p_fp", "p_visit")
  miss <- setdiff(need, names(sp))
  if (length(miss)) stop("species table lacks columns: ", paste(miss, collapse = ", "))
  sp$visitor_only <- as.logical(sp$visitor_only)
  sp$confusable <- if ("confusable" %in% names(sp)) as.logical(sp$confusable) else FALSE
  if (!"frag_affected" %in% names(sp)) {          # derive from the betas if the column is missing
    bc <- setdiff(grep("^beta_", names(sp), value = TRUE), "beta0")
    sp$frag_affected <- ifelse(sp$visitor_only, NA, rowSums(abs(sp[bc]), na.rm = TRUE) > 0)
  }
  sp$frag_affected <- as.logical(sp$frag_affected)
  if (!"beta0" %in% names(sp)) sp$beta0 <- NA_real_
  sp
}

# ===========================================================================
# Latent occupancy / breeding state ####
# ===========================================================================
# Latent occupancy / breeding state: z_sc ~ Bernoulli(psi_sc),
#   logit psi_sc = beta0_s + sum_k beta_sk x_ck + field_s(c)
# beta0_s is solved (if NA) so that the mean psi over non-water cells equals the species'
# target prevalence. Visitor-only species (if any) never occupy (psi = 0, z = 0).
#
# Spatial field: species with a non-NA field_scale_m get a Gaussian random field with
# exponential correlation exp(-d / field_scale_m) and SD cfg$occupancy$field_sd (logit scale).
# If cfg$occupancy$field_standardise is TRUE the field is rescaled so that its SD over land
# cells is exactly field_sd (and its land mean is 0). This matters for large scales: on a
# 50 km grid a 10 km field would otherwise express only part of its nominal variance, and
# its mean would be absorbed by beta0.
#
# Returns list(species, psi, z, field):
#   psi, z : arrays [nrow, ncol, n_species]
#   field  : list of n_species matrices (field SD already applied), NULL where no field

sim_build_occupancy <- function(cfg, land, sp) {
  nr <- land$nrow
  nc <- land$ncol
  n <- nrow(sp)
  ok <- !as.vector(land$water)
  bcols <- paste0("beta_", cfg$cov_names)
  if (!"beta0" %in% names(sp)) sp$beta0 <- NA_real_
  fsd <- cfg$occupancy$field_sd
  fstd <- isTRUE(cfg$occupancy$field_standardise)

  psi <- array(0, dim = c(nr, nc, n))
  z <- array(0L, dim = c(nr, nc, n))
  field <- vector("list", n)
  setups <- list()   # one FFT setup per distinct scale, reused across species

  for (s in seq_len(n)) {
    if (isTRUE(sp$visitor_only[s])) next

    beta <- as.numeric(unlist(sp[s, bcols]))
    eta <- as.vector(land$X %*% beta)

    scale_s <- sp$field_scale_m[s]
    if (!is.na(scale_s) && fsd > 0) {
      key <- as.character(scale_s)
      if (is.null(setups[[key]])) {
        setups[[key]] <- sim_grf_setup(nr, nc, cfg$cell_m, scale_s)
        if (setups[[key]]$neg_frac > 0.01)
          warning("field scale ", scale_s, " m: ", round(100 * setups[[key]]$neg_frac, 1),
                  "% of the spectrum was clipped; the field is only approximately exponential")
      }
      g <- sim_with_seed(sim_seed(cfg$seed, 11L, s), sim_grf_draw(setups[[key]]))
      if (fstd) g <- sim_standardise(g, !land$water)
      field[[s]] <- fsd * g
      eta <- eta + as.vector(field[[s]])
    }

    if (is.na(sp$beta0[s])) {
      target <- sp$prevalence_target[s]
      f <- function(b) mean(plogis(b + eta[ok])) - target
      sp$beta0[s] <- uniroot(f, interval = c(-10, 5), extendInt = "yes")$root
    }

    p <- plogis(sp$beta0[s] + eta)
    p[!ok] <- 0
    u <- sim_with_seed(sim_seed(cfg$seed, 12L, s), runif(nr * nc))
    psi[, , s] <- matrix(p, nr, nc)
    z[, , s] <- matrix(as.integer(u < p), nr, nc)
  }
  list(species = sp, psi = psi, z = z, field = field)
}

# ===========================================================================
# Phenology: Acoustic activity curve
# ===========================================================================
# Acoustic activity curve a(d) in [0, 1]: 0 outside [onset, cessation], 1 at the peak.
# Shape: Beta density on u = (d - onset) / (cessation - onset), with mode at
# m = (peak - onset) / (cessation - onset) and concentration kappa
# (Beta(1 + kappa m, 1 + kappa (1 - m))), scaled to a maximum of 1.
# Vectorised: d, onset, peak, cessation, kappa may be vectors (recycled).

sim_activity <- function(d, onset, peak, cessation, kappa = 4) {
  n <- max(length(d), length(onset), length(peak), length(cessation), length(kappa))
  d <- rep_len(d, n)
  onset <- rep_len(onset, n)
  peak <- rep_len(peak, n)
  cessation <- rep_len(cessation, n)
  kappa <- rep_len(kappa, n)

  w <- cessation - onset
  u <- (d - onset) / w
  m <- (peak - onset) / w
  out <- numeric(n)
  ins <- !is.na(u) & u > 0 & u < 1 & !is.na(m) & m > 0 & m < 1
  if (any(ins)) {
    ui <- u[ins]
    mi <- m[ins]
    ki <- kappa[ins]
    out[ins] <- exp(ki * mi * (log(ui) - log(mi)) + ki * (1 - mi) * (log(1 - ui) - log(1 - mi)))
  }
  out
}

# Matrix [n_species, length(doys)] of species activity (no cell shift), for plots and checks.
sim_activity_curves <- function(sp, doys) {
  t(vapply(seq_len(nrow(sp)), function(s)
    sim_activity(doys, sp$onset[s], sp$peak[s], sp$cessation[s], sp$shape_kappa[s]),
    numeric(length(doys))))
}

# ===========================================================================
# Observation model for one 24 h recording day at cell (i, j) on day "doy" ####
# ===========================================================================
# The recorder takes fixed 1-minute snippets (cfg$snippets): one every interval_min minutes,
# starting at minute 0 of the day (e.g. every 10 min = 144 snippets per day). Calling rates
# are flat over the 24 h, so every snippet has the same chance of a call.
#
# Three independent sources per species, drawn snippet by snippet:
#   1. true calls  (z = 1): calls per snippet ~ Poisson(lambda x activity(doy - cell shift) x length / 1440),
#                           each call scores ~ U(0.8, 1)
#   2. visitor     (z = 0): the species visits today with probability p_visit x activity(doy); while it
#                           visits it calls at lambda x activity(doy) x visit_frac per day (Poisson per
#                           snippet), each call scores ~ U(0.6, 1)
#   3. false pos.  (always): each snippet has a false positive with probability p_fp (per snippet),
#                           score ~ U(0.3, 1)
# The maximum score of N calls on [lo, hi] is drawn directly: lo + (hi - lo) U^(1/N).
#
# Outputs per species:
#   max_score    highest score over all snippets (0 = no detection)
#   n_snip_det   number of snippets with at least one detection of any source, whatever the score
#
# Reproducible: all random numbers come from one seed (seed, i, j, doy), so results do not
# depend on what was queried before.
#
# sim_cell_day() returns a data.frame with one row per species:
#   species_id, max_score, n_snip_det, n_true (true calls in the snippets), n_true_snips,
#   n_vis_snips, n_fp_snips (snippets with each source), visitor (visitor called in a snippet),
#   fp (any false positive), source of the max ("none", "true", "visitor", "fp"),
#   z (true state), act (cell activity, 0-1).
# simulate() exposes only species_id, max_score and n_snip_det; the rest is truth.

sim_n_snippets <- function(cfg) {
  sn <- cfg$snippets
  if (sn$length_min > sn$interval_min) stop("snippets$length_min must not exceed interval_min")
  floor(1440 / sn$interval_min)
}

# Start minute (0-1439) of every snippet of the day. Fixed schedule.
sim_snippet_starts <- function(cfg) (seq_len(sim_n_snippets(cfg)) - 1) * cfg$snippets$interval_min

sim_cell_day <- function(world, i, j, doy) {
  cfg <- world$cfg
  sp <- world$species
  n <- nrow(sp)
  ob <- cfg$obs
  S <- sim_n_snippets(cfg)
  per_snip <- cfg$snippets$length_min / 1440   # share of the day covered by one snippet

  # each species' timing (phenology) shift in this cell, drawn from Normal(0, cell_shift_sd)
  delta <- sim_with_seed(sim_seed(cfg$seed, 13L, i, j),
                         rnorm(n, 0, cfg$pheno$cell_shift_sd))
  # true state of each species in the cell
  z <- as.integer(world$z[i, j, ])
  # activity (0 to 1) of each species today, in this cell (shifted curve) and in general (visitors)
  act_cell <- sim_activity(doy - delta, sp$onset, sp$peak, sp$cessation, sp$shape_kappa)
  act_sp <- sim_activity(doy, sp$onset, sp$peak, sp$cessation, sp$shape_kappa)

  # mean true calls per snippet (only where the species occupies the cell)
  mu_true <- ifelse(z == 1L, sp$lambda * act_cell, 0) * per_snip
  mu_true[is.na(mu_true)] <- 0

  draw <- sim_with_seed(sim_seed(cfg$seed, 20L, i, j, doy), {
    # visitor today? only where the species does not occupy the cell
    visitor_day <- (z == 0L) & (runif(n) < sp$p_visit * act_sp)
    mu_vis <- ifelse(visitor_day, sp$lambda * act_sp * cfg$visitor$visit_frac, 0) * per_snip
    mu_vis[is.na(mu_vis)] <- 0
    list(true = matrix(rpois(n * S, rep(mu_true, S)), n, S),     # calls per snippet, species x snippet
         vis  = matrix(rpois(n * S, rep(mu_vis, S)), n, S),
         fp   = matrix(runif(n * S) < rep(sp$p_fp, S), n, S),     # false positive in the snippet?
         u    = matrix(runif(3 * n), n, 3))                       # uniforms for the three max scores
  })

  n_true <- rowSums(draw$true)
  n_vis <- rowSums(draw$vis)
  n_fp <- rowSums(draw$fp)

  max_true <- ifelse(n_true > 0,
                     ob$true_score[1] + diff(ob$true_score) * draw$u[, 1]^(1 / pmax(n_true, 1)), 0)
  max_vis <- ifelse(n_vis > 0,
                    ob$visitor_score[1] + diff(ob$visitor_score) * draw$u[, 2]^(1 / pmax(n_vis, 1)), 0)
  max_fp <- ifelse(n_fp > 0,
                   ob$fp_score[1] + diff(ob$fp_score) * draw$u[, 3]^(1 / pmax(n_fp, 1)), 0)

  M <- cbind(max_true, max_vis, max_fp)
  mx <- pmax(max_true, max_vis, max_fp)
  which_max <- max.col(M, ties.method = "first")
  source <- c("none", "true", "visitor", "fp")[1L + ifelse(mx > 0, which_max, 0L)]

  data.frame(species_id = sp$species_id, max_score = mx,
             n_snip_det = rowSums(draw$true > 0 | draw$vis > 0 | draw$fp),
             n_true = n_true,
             n_true_snips = rowSums(draw$true > 0), n_vis_snips = rowSums(draw$vis > 0),
             n_fp_snips = n_fp,
             visitor = n_vis > 0, fp = n_fp > 0, source = source, z = z, act = act_cell,
             stringsAsFactors = FALSE)
}
