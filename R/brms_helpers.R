# Shared data preparation, priors and helpers for the brms models
# (analysis/01_models.R, 02_initial_composition.R, 03_disturbance.R)

z <- \(x) as.numeric(scale(x))

# Standardize a composition-level predictor within richness levels:
# (x - mean within div) / pooled within-div SD. Slopes then measure the effect
# of compositional differences at a fixed richness.
z_within <- function(x, g) {
  dev <- x - ave(x, g)
  dev / sqrt(sum(dev^2) / (length(dev) - length(unique(g))))
}

SEED <- 20261015
ITER <- 2000
ITER_MV <- 4000
CONTROL <- list(adapt_delta = 0.95)

responses <- c(prod = "y_prod", resp = "y_resp", dc = "y_dc")
# prod on the raw scale (log(prod) is left-skewed), log(resp), and dc with
# Student-t errors (heavy tails), as indicated by posterior predictive checks
families <- list(y_prod = gaussian(), y_resp = gaussian(), y_dc = student())
# T98: residual SD of prod and resp varies among sampling periods (grouped PPC)
sigma98 <- list(y_prod = "ftime", y_resp = "ftime")

prep_data <- function(bef) {
  bef |>
    dplyr::mutate(
      y_prod = z(prod),
      y_resp = z(log(resp)),
      y_dc   = z(dc),
      div    = factor(div, levels = c(1, 9, 15), labels = c("S3", "S9", "S15")),
      ftime  = factor(time),
      comp   = factor(comp),
      tank   = factor(tank)
    )
}

# ---- Priors (standardized scale) ----
priors_uni <- c(
  prior(normal(0, 1), class = b),
  prior(exponential(2), class = sd),
  prior(exponential(1), class = sigma),
  prior(normal(0, 0.5), class = ar)
)
prior_nu <- prior(gamma(2, 0.1), class = nu)

# x = NULL gives the baseline model without a phylogenetic predictor;
# `extra` adds further fixed-effect terms
rhs <- function(x = NULL, ar = TRUE, extra = NULL) {
  paste0(
    "0 + div", if (!is.null(x)) paste0(" + div:", x),
    " + ftime + dist + dist:post", if (!is.null(extra)) paste0(" + ", extra),
    " + (1 | comp) + (1 | tank)",
    if (ar) " + ar(time = time, gr = tank, p = 1)" else ""
  )
}

# Distributional model for the residual SD (log scale; standardized responses)
priors_sigma <- c(
  prior(normal(0, 1), class = Intercept, dpar = sigma),
  prior(normal(0, 1), class = b, dpar = sigma)
)

# `rhs_override` replaces the default right-hand side (used in step 3);
# `sigma_rhs` (e.g. "ftime") models the residual SD
fit_uni <- function(resp, x = NULL, data, ar = TRUE, extra = NULL, name = NULL,
                    rhs_override = NULL, sigma_rhs = NULL) {
  if (is.null(name)) name <- paste(c("uni", resp, x %||% "base", if (ar) "ar1" else "noar",
                                     if (!is.null(sigma_rhs)) "sig"), collapse = "_")
  fml <- paste(resp, "~", rhs_override %||% rhs(x, ar, extra))
  has_ar <- ar && grepl("ar(", fml, fixed = TRUE)
  pri <- if (has_ar) priors_uni else priors_uni[priors_uni$class != "ar", ]
  if (!is.null(sigma_rhs)) pri <- c(pri[pri$class != "sigma", ], priors_sigma)
  args <- list(
    formula = if (is.null(sigma_rhs)) bf(fml) else bf(as.formula(fml), as.formula(paste("sigma ~", sigma_rhs))),
    data = data, family = families[[resp]],
    prior = c(pri, if (families[[resp]]$family == "student") prior_nu),
    iter = ITER, chains = 4, seed = SEED, control = CONTROL,
    file = here::here("output", "models", name), file_refit = "on_change"
  )
  fit <- do.call(brm, args)
  # Refit with twice the iterations and adapt_delta = 0.99 when convergence
  # criteria are not met (R-hat > 1.01 or bulk/tail ESS < 400)
  dg <- diagnose(fit, name)
  if (dg$max_rhat > 1.01 || min(dg$min_ess_bulk, dg$min_ess_tail) < 400) {
    message("Refitting ", name, " with more iterations")
    args$iter <- 2 * ITER
    args$control <- list(adapt_delta = 0.99)
    args$file <- here::here("output", "models", paste0(name, "_long"))
    fit <- do.call(brm, args)
  }
  fit
}

# ---- Diagnostics ----
diagnose <- function(fit, name) {
  s <- posterior::summarise_draws(posterior::as_draws_df(fit), "rhat", "ess_bulk", "ess_tail")
  np <- nuts_params(fit)
  tibble::tibble(
    model = name,
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(s$ess_tail, na.rm = TRUE),
    divergences = sum(np$Value[np$Parameter == "divergent__"]),
    max_treedepth_hits = sum(np$Value[np$Parameter == "treedepth__"] >= 10)
  )
}

# Slopes of predictor x per richness level, plus their average
slope_draws <- function(fit, x, resp_prefix = "", levels = c("S3", "S9", "S15")) {
  dr <- posterior::as_draws_df(fit)
  cols <- paste0("b_", resp_prefix, "div", levels, ":", x)
  out <- purrr::map(cols, \(cl) dr[[cl]]) |> purrr::set_names(levels)
  out$average <- rowMeans(do.call(cbind, out))
  out
}

# Posterior summary with 95% and Bonferroni-type 98.33% intervals
summarise_draws_ci <- function(v) {
  tibble::tibble(
    mean = mean(v),
    q2.5 = quantile(v, 0.025), q97.5 = quantile(v, 0.975),
    q0.83 = quantile(v, 0.05 / 6), q99.17 = quantile(v, 1 - 0.05 / 6),
    p_gt0 = mean(v > 0)
  )
}

summarise_slopes <- function(sl, model, response, x) {
  purrr::imap_dfr(sl, \(v, k) dplyr::bind_cols(
    tibble::tibble(model = model, response = response, predictor = x, div = k),
    summarise_draws_ci(v)
  ))
}

save_ppc <- function(fit, name, resp = NULL) {
  p <- pp_check(fit, ndraws = 100, resp = resp) + ggplot2::ggtitle(name)
  ggplot2::ggsave(here::here("figures", "ppc", paste0(name, ".png")), p,
                  width = 5, height = 3.5, dpi = 150)
}
