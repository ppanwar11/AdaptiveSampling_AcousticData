# ===========================================================================
# Functions to simulated data ####
# ===========================================================================
# Public interface.
#   sim_init(cfg)              build the world (landscape, species, occupancy), keep it as default
#   simulate(i, j, doy)        one recording day: data.frame(species_id, max_score, n_snip_det)
#   sim_collect(design)        many recordings; design has columns i, j, doy
#   sim_truth(), sim_truth_day(i, j, doy), sim_truth_maps()   TRUTH, not for the sampler
#   sim_plot_truth(species_id)
# i = row (1 = north), j = column (1 = west). doy in [doy_min, doy_max].
# NB: simulate() masks stats::simulate; call stats::simulate() if you need that one.

sim_init <- function(cfg = sim_config(), species = NULL, set_default = TRUE) {
  land <- sim_load_covariates(cfg)
  sp <- if (is.null(species)) sim_make_species(cfg, land) else species
  occ <- sim_build_occupancy(cfg, land, sp)
  world <- structure(list(cfg = cfg, land = land, species = occ$species,
                          psi = occ$psi, z = occ$z, field = occ$field),
                     class = "sim_world")
  if (set_default) assign("world", world, envir = .sim_env)
  invisible(world)
}

sim_get_world <- function(world = NULL) {
  if (!is.null(world)) return(world)
  w <- get0("world", envir = .sim_env, inherits = FALSE)
  if (is.null(w)) stop("no world yet: run sim_init() first")
  w
}

sim_check_query <- function(w, i, j, doy) {
  if (length(i) != 1 || length(j) != 1 || length(doy) != 1) stop("i, j and doy must be single numbers")
  if (anyNA(c(i, j, doy))) stop("i, j and doy must not be NA")
  if (i != round(i) || j != round(j) || doy != round(doy)) stop("i, j and doy must be whole numbers")
  if (i < 1 || i > w$land$nrow || j < 1 || j > w$land$ncol) stop("(i, j) outside the grid")
  if (doy < w$cfg$doy_min || doy > w$cfg$doy_max)
    stop("doy must be between ", w$cfg$doy_min, " and ", w$cfg$doy_max)
  if (w$land$water[i, j]) stop("cell (", i, ", ", j, ") is excluded (water)")
  invisible(TRUE)
}

# One recording day (fixed 1-minute snippets, see cfg$snippets). Returns species_id, max_score
# (highest score over the day, 0 = no detection) and n_snip_det (number of snippets with any
# detection, whatever the score).
simulate <- function(i, j, doy, world = NULL) {
  w <- sim_get_world(world)
  sim_check_query(w, i, j, doy)
  sim_cell_day(w, i, j, doy)[, c("species_id", "max_score", "n_snip_det")]
}

# Many recordings. design: data.frame with columns i, j, doy.
sim_collect <- function(design, world = NULL) {
  w <- sim_get_world(world)
  out <- vector("list", nrow(design))
  for (r in seq_len(nrow(design))) {
    d <- simulate(design$i[r], design$j[r], design$doy[r], world = w)
    out[[r]] <- cbind(design[rep(r, nrow(d)), c("i", "j", "doy"), drop = FALSE], d)
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

# ---- truth accessors -------------------------------------------------------
sim_truth <- function(world = NULL) {
  w <- sim_get_world(world)
  list(cfg = w$cfg, species = w$species, psi = w$psi, z = w$z, field = w$field,
       water = w$land$water, covariates = w$land$x)
}

# Full latent detail for one recording: true calls, visitor, FP, source of the max.
sim_truth_day <- function(i, j, doy, world = NULL) {
  w <- sim_get_world(world)
  sim_check_query(w, i, j, doy)
  sim_cell_day(w, i, j, doy)
}

# Maps: richness (sum of z), mean psi; for one species also its psi, z and field.
sim_truth_maps <- function(world = NULL, species_id = NULL) {
  w <- sim_get_world(world)
  res <- list(richness = apply(w$z, c(1, 2), sum),
              mean_psi = apply(w$psi, c(1, 2), mean))
  res$richness[w$land$water] <- NA
  res$mean_psi[w$land$water] <- NA
  if (!is.null(species_id)) {
    s <- match(species_id, w$species$species_id)
    if (is.na(s)) stop("unknown species_id")
    res$psi <- w$psi[, , s]
    res$z <- w$z[, , s]
    res$field <- w$field[[s]]
    res$psi[w$land$water] <- NA
    res$z[w$land$water] <- NA
  }
  res
}

# Draw a matrix with north up and west left.
sim_image <- function(m, main = "", col = hcl.colors(50, "YlGnBu", rev = TRUE), ...) {
  image(x = seq_len(ncol(m)), y = seq_len(nrow(m)), z = t(m), ylim = c(nrow(m), 1),
        asp = 1, col = col, main = main, xlab = "j (west to east)", ylab = "i (north to south)", ...)
}

sim_plot_truth <- function(species_id, world = NULL) {
  w <- sim_get_world(world)
  m <- sim_truth_maps(w, species_id)
  op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  on.exit(par(op))
  sim_image(m$psi, main = paste(species_id, "psi"))
  sim_image(m$z, main = paste(species_id, "occupied (z)"), col = c("grey90", "darkgreen"))
  sim_image(m$richness, main = "richness (true)")
  invisible(m)
}

# ===========================================================================
# Checks that the simulated world behaves as designed. 
# ===========================================================================
# Prints tables, optional plots, returns everything invisibly as a list.
#   v <- sim_validate()                       # uses the default world
#   v <- sim_validate(n_cells = 300, n_days = 30, plots = FALSE)

sim_validate <- function(world = NULL, n_cells = 200, n_days = 20, plots = TRUE, seed = 1) {
  w <- sim_get_world(world)
  sp <- w$species
  cfg <- w$cfg
  nr <- w$land$nrow
  nc <- w$land$ncol
  ok <- !w$land$water
  idx <- which(ok)
  k <- length(cfg$cov_names)
  bcols <- paste0("beta_", cfg$cov_names)
  occ_idx <- which(!sp$visitor_only)
  out <- list()

  # 1. prevalence ------------------------------------------------------------
  prev <- apply(w$z, 3, function(m) mean(m[ok]))
  out$prevalence <- data.frame(species_id = sp$species_id, strategy = sp$strategy,
                               profile = sp$profile, target = sp$prevalence_target,
                               realised = prev)
  cat("\n[1] Prevalence: realised vs target (occupying species), correlation =",
      round(cor(prev[occ_idx], sp$prevalence_target[occ_idx]), 3), "\n")
  if (any(sp$visitor_only)) {
    cat("    visitor-only species max prevalence:", max(prev[sp$visitor_only]), "(should be 0)\n")
  } else {
    cat("    no visitor-only species in this world\n")
  }

  sim_with_seed(seed, {
    # 2. effect recovery -------------------------------------------------------
    sub <- sample(idx, min(20000, length(idx)))
    Xs <- as.data.frame(w$land$X[sub, , drop = FALSE])
    names(Xs) <- cfg$cov_names
    est <- t(vapply(occ_idx, function(s) {
      y <- w$z[, , s][sub]
      if (sum(y) < 10 || sum(1 - y) < 10) return(rep(NA_real_, k))
      dat <- Xs
      dat$y <- y
      fit <- suppressWarnings(stats::glm(y ~ ., data = dat, family = stats::binomial))
      unname(coef(fit)[cfg$cov_names])
    }, numeric(k)))
    tru <- as.matrix(sp[occ_idx, bcols, drop = FALSE])
    aff <- if ("frag_affected" %in% names(sp)) sp$frag_affected[occ_idx] %in% TRUE else rep(TRUE, length(occ_idx))
    rec <- do.call(rbind, lapply(seq_len(k), function(q) {
      ok_q <- !is.na(est[, q]) & aff                         # affected species only
      data.frame(covariate = cfg$cov_names[q],
                 config_mean = cfg$occupancy$beta_mean[[cfg$cov_names[q]]],
                 config_sd = cfg$occupancy$beta_sd[[cfg$cov_names[q]]],
                 drawn_mean = mean(tru[aff, q]), drawn_sd = sd(tru[aff, q]),
                 cor_est_true = cor(est[ok_q, q], tru[ok_q, q]),
                 slope_est_on_true = unname(coef(lm(est[ok_q, q] ~ tru[ok_q, q]))[2]))
    }))
    out$effect_recovery <- rec
    cat("\n[2] Covariate effects of the species affected by fragmentation (", sum(aff), " of ",
        length(aff), " occupying species): drawn betas vs config, and logistic-regression recovery\n", sep = "")
    cat("    (slope < 1 is expected: spatial fields add noise and attenuate estimates)\n")
    print(rec, digits = 3, row.names = FALSE)

    # species not affected by fragmentation: true betas are 0; what does the regression estimate?
    un <- !aff
    if (any(un)) {
      out$unaffected_estimates <- data.frame(
        covariate = cfg$cov_names,
        n_species = sum(un & !is.na(est[, 1])),
        mean_abs_estimate = vapply(seq_len(k), function(q) mean(abs(est[un, q]), na.rm = TRUE), numeric(1)),
        max_abs_estimate = vapply(seq_len(k), function(q) max(abs(est[un, q]), na.rm = TRUE), numeric(1)),
        share_abs_gt_0.2 = vapply(seq_len(k), function(q) mean(abs(est[un, q]) > 0.2, na.rm = TRUE), numeric(1)))
      cat("\n[2b] Species NOT affected by fragmentation (true betas = 0): estimated effects from the same regression\n")
      cat("     (non-zero estimates come from the spatial fields and covariate correlation)\n")
      print(out$unaffected_estimates, digits = 3, row.names = FALSE)
    }

    # 3. spatial autocorrelation of the fields ----------------------------------
    with_field <- which(!is.na(sp$field_scale_m))
    out$field_scales <- table(field_scale_m = sp$field_scale_m, useNA = "ifany")
    cat("\n[3] Species per field e-folding distance (NA = no spatial field)\n")
    print(out$field_scales)
    if (length(with_field)) {
      okm <- !w$land$water
      lags <- c(1, 2, 5, 10, 20, 50)
      lags <- lags[lags < nc]
      show <- with_field[!duplicated(sp$field_scale_m[with_field])]   # one species per distinct scale
      out$autocorrelation <- do.call(rbind, lapply(show, function(s) {
        f <- w$field[[s]] / cfg$occupancy$field_sd
        emp <- vapply(lags, function(l) {
          a <- f[, 1:(nc - l)]
          b <- f[, (1 + l):nc]
          both <- okm[, 1:(nc - l)] & okm[, (1 + l):nc]       # land pairs only
          cor(a[both], b[both])
        }, numeric(1))
        data.frame(species_id = sp$species_id[s], scale_m = sp$field_scale_m[s],
                   field_sd_land = sd(w$field[[s]][okm]),
                   lag_m = lags * cfg$cell_m, empirical = emp,
                   expected = exp(-lags * cfg$cell_m / sp$field_scale_m[s]))
      }))
      cat("\n    Field autocorrelation over land cell pairs, empirical vs exp(-d / scale), one species per scale\n")
      cat("    (a single field on a finite grid gives empirical values below expected at large scales)\n")
      print(out$autocorrelation, digits = 2, row.names = FALSE)
    }
    cat("    species with a spatial field:", length(with_field), "of", nrow(sp), "\n")

    # 4. phenology -----------------------------------------------------------------
    out$phenology <- stats::aggregate(cbind(onset, peak, cessation) ~ strategy, data = sp, FUN = mean)
    cat("\n[4] Mean onset / peak / cessation (DOY) by strategy\n")
    print(out$phenology, digits = 3, row.names = FALSE)

    # 5. simulated recordings -----------------------------------------------------------
    n_cells <- min(n_cells, length(idx))
    cells <- idx[sample.int(length(idx), n_cells)]
    ci <- ((cells - 1) %% nr) + 1
    cj <- ((cells - 1) %/% nr) + 1
    days <- cfg$doy_min:cfg$doy_max
    rows <- vector("list", n_cells * n_days)
    r <- 0
    for (cc in seq_len(n_cells)) {
      for (d in sample(days, n_days)) {
        r <- r + 1
        x <- sim_cell_day(w, ci[cc], cj[cc], d)
        x$i <- ci[cc]
        x$j <- cj[cc]
        x$doy <- d
        rows[[r]] <- x
      }
    }
    cd <- do.call(rbind, rows)
    m <- match(cd$species_id, sp$species_id)
    cd$profile <- sp$profile[m]
    cd$strategy <- sp$strategy[m]
    cd$confusable <- sp$confusable[m]
    cd$group <- ifelse(cd$confusable, paste0(cd$profile, "+conf"), cd$profile)   # profile, confusable marked
    cd$detected <- cd$max_score > 0
    out$recordings <- cd

    # 5a. detection of occupying species near their peak
    pk <- cd[cd$z == 1 & cd$act >= 0.8, ]
    if (nrow(pk)) {
      out$detect_peak <- data.frame(
        group = sort(unique(pk$group)),
        p_true_call = as.vector(tapply(pk$n_true > 0, pk$group, mean)),
        mean_max_score = as.vector(tapply(ifelse(pk$detected, pk$max_score, NA), pk$group,
                                          mean, na.rm = TRUE)),
        mean_snip_det = as.vector(tapply(pk$n_snip_det, pk$group, mean)),
        n = as.vector(table(pk$group)))
      cat("\n[5a] Occupied cells at activity >= 0.8 (", sim_n_snippets(cfg), " snippets per day): ",
          "probability of a true call in the day's snippets, mean max score, mean number of snippets with a detection\n", sep = "")
      print(out$detect_peak, digits = 3, row.names = FALSE)
    }

    # 5b. false positives and visitors where the species does not occupy
    ab <- cd[cd$z == 0, ]
    out$absent <- data.frame(
      group = sort(unique(ab$group)),
      p_fp = as.vector(tapply(ab$fp, ab$group, mean)),
      p_fp_ge_0.8 = as.vector(tapply(ab$fp & ab$max_score >= 0.8 & ab$source == "fp", ab$group, mean)),
      p_visitor = as.vector(tapply(ab$visitor, ab$group, mean)),
      p_any_detection = as.vector(tapply(ab$detected, ab$group, mean)))
    cat("\n[5b] Cell-days where the species does not occupy: FP, visitor and any detection rates\n")
    print(out$absent, digits = 3, row.names = FALSE)

    # 5c. what produces scores >= 0.8
    hi <- cd[cd$max_score >= 0.8, ]
    out$source_ge_0.8 <- table(group = hi$group, source = hi$source)
    cat("\n[5c] Source of the max score when it is >= 0.8\n")
    print(out$source_ge_0.8)

    # 5d. calibration: is the score a probability?
    det <- cd[cd$detected, ]
    det$bin <- cut(det$max_score, breaks = seq(0.3, 1, by = 0.1), include.lowest = TRUE)
    det$correct <- det$source %in% c("true", "visitor")
    mids <- seq(0.35, 0.95, by = 0.1)
    out$calibration <- data.frame(bin = levels(det$bin), score_mid = mids,
                                  n = as.vector(table(det$bin)),
                                  p_correct = as.vector(tapply(det$correct, det$bin, mean)))
    cat("\n[5d] Calibration: share of max scores that are real bird calls (true or visitor)\n")
    cat("     If the score were an exact probability, p_correct would equal score_mid.\n")
    print(out$calibration, digits = 3, row.names = FALSE)

    # 5d2. snippet counts agree with the scores
    out$snip_consistent <- all((cd$max_score > 0) == (cd$n_snip_det > 0)) &&
      all(cd$n_snip_det <= sim_n_snippets(cfg))
    cat("\n[5d2] n_snip_det > 0 exactly when max_score > 0, and never above the number of snippets:",
        out$snip_consistent, "\n")

    # 5e. visitor-only species
    vo <- cd[cd$profile == "visitor_only", ]
    if (nrow(vo)) {
      out$visitor_only <- data.frame(
        species_id = sort(unique(vo$species_id)),
        p_detected = as.vector(tapply(vo$detected, vo$species_id, mean)),
        p_visitor = as.vector(tapply(vo$visitor, vo$species_id, mean)))
      cat("\n[5e] Visitor-only species: share of recordings with any detection\n")
      print(out$visitor_only, digits = 3, row.names = FALSE)
    }
  })

  # 6. reproducibility ----------------------------------------------------------------
  ci1 <- ((idx[1] - 1) %% nr) + 1
  cj1 <- ((idx[1] - 1) %/% nr) + 1
  set.seed(99)
  rng0 <- .Random.seed
  a <- simulate(ci1, cj1, cfg$doy_min + 30, world = w)
  simulate(ci1, cj1, cfg$doy_min + 31, world = w)
  b <- simulate(ci1, cj1, cfg$doy_min + 30, world = w)
  out$reproducible <- identical(a, b)
  out$rng_untouched <- identical(rng0, .Random.seed)
  cat("\n[6] Same query gives identical output after other queries:", out$reproducible,
      "| caller's RNG state untouched:", out$rng_untouched, "\n")

  # plots -------------------------------------------------------------------------------
  if (plots) {
    op <- par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))
    on.exit(par(op), add = TRUE)
    doys <- 60:260
    curves <- sim_activity_curves(sp, doys)
    cols <- c(resident = "#1b9e77", short = "#d95f02", long = "#7570b3")[sp$strategy]
    plot(NA, xlim = range(doys), ylim = c(0, 1), xlab = "day of year", ylab = "activity",
         main = "Activity curves by strategy")
    for (s in seq_len(nrow(sp))) lines(doys, curves[s, ], col = cols[s])
    abline(v = c(cfg$doy_min, cfg$doy_max), lty = 2, col = "grey50")
    legend("topright", legend = c("resident", "short", "long"),
           col = c("#1b9e77", "#d95f02", "#7570b3"), lty = 1, bty = "n", cex = 0.8)
    sim_image(sim_truth_maps(w)$richness, main = "True species richness")
    cd0 <- out$recordings
    hist(cd0$max_score[cd0$detected], breaks = seq(0.3, 1, by = 0.025), col = "grey70",
         main = "Max scores (detections only)", xlab = "max score")
    # plot(out$calibration$score_mid, out$calibration$p_correct, type = "b", pch = 16,
    #      xlim = c(0.3, 1), ylim = c(0, 1), xlab = "max score", ylab = "share real calls",
    #      main = "Calibration")
    # abline(0, 1, lty = 2, col = "grey50")
  }
  invisible(out)
}
