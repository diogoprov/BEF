# Bayesian check of the design simulation (analysis/18_design_simulation.R).
#
# The REML Wald test in script 18 is anti-conservative for D2/D3 at
# n_comp = 24, plausibly because the variance components are poorly estimated
# near zero. Here the same generative model is fitted with brms (posterior
# propagates variance-component uncertainty) for D2 and D4 (first simulated
# tree), n_comp = 24, theta = 0, scenario phy+id. Reported: rate at which the
# 95% CrI of theta excludes zero (target 0.05), and MCMC diagnostics.
#
# Inputs: output/models/design_sim_designs.rds (written by script 18)
# Outputs: output/tables/design_sim_brms.csv, design_sim_brms_diagnostics.csv;
#          output/models/design_sim_brms_reps.rds (git-ignored)

suppressPackageStartupMessages({library(tidyverse); library(brms); library(posterior); library(here)})
options(mc.cores = 4, brms.backend = "cmdstanr")
source(here("R", "brms_helpers.R"))

N_REP <- 100
S_PHY <- 0.15; S_ID <- 0.15; S_COMP <- 0.10; S_TANK <- 0.30; N_TANK <- 4
set.seed(SEED)

designs <- readRDS(here("output", "models", "design_sim_designs.rds"))[["24"]]
designs <- list(D2 = designs$D2[[1]], D4 = designs$D4[[1]])

rmvn <- function(K) {
  e <- eigen(K, symmetric = TRUE)
  drop(e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(nrow(K))))
}

# Projected kernels are rank-deficient; jitter as in script 13
jitter_k <- function(K, ids) {
  K <- K + diag(1e-4, nrow(K)); dimnames(K) <- list(ids, ids); K
}

priors_sim <- c(
  prior(normal(0, 1), class = b),
  prior(normal(0, 1), class = Intercept),
  prior(exponential(2), class = sd),
  prior(exponential(2), class = sigma)
)

sim_data <- function(des) {
  n <- nrow(des$P)
  ids <- as.character(seq_len(n))
  zz <- (des$mpd - mean(des$mpd)) / des$ref_sd
  y <- S_PHY * rmvn(des$K_phy) + S_ID * rmvn(des$K_id) +
    rnorm(n, 0, sqrt(S_COMP^2 + S_TANK^2 / N_TANK))
  tibble(y = y, zmpd = zz, c_phy = ids, c_id = ids)
}

run_design <- function(nm) {
  des <- designs[[nm]]
  ids <- as.character(seq_len(nrow(des$P)))
  d2 <- list(K_phy = jitter_k(des$K_phy, ids), K_id = jitter_k(des$K_id, ids))
  base <- brm(y ~ 1 + zmpd + (1 | gr(c_phy, cov = K_phy)) + (1 | gr(c_id, cov = K_id)),
              data = sim_data(des), data2 = d2, prior = priors_sim,
              chains = 4, iter = ITER, seed = SEED, control = list(adapt_delta = 0.99),
              refresh = 0, silent = 2)
  map_dfr(seq_len(N_REP), \(r) {
    fit <- update(base, newdata = sim_data(des), data2 = d2, recompile = FALSE,
                  seed = SEED + r, refresh = 0, silent = 2)
    dr <- as_draws_df(fit)
    np <- nuts_params(fit)
    pars <- c("b_zmpd", "sd_c_phy__Intercept", "sd_c_id__Intercept", "sigma")
    sm <- summarise_draws(subset_draws(dr, variable = pars), "rhat", "ess_bulk", "ess_tail")
    tibble(
      design = nm, rep = r,
      est = median(dr$b_zmpd),
      lo = quantile(dr$b_zmpd, 0.025), hi = quantile(dr$b_zmpd, 0.975),
      divergences = sum(np$Value[np$Parameter == "divergent__"]),
      diag = list(sm)
    )
  })
}

res <- map_dfr(names(designs), run_design)

saveRDS(res, here("output", "models", "design_sim_brms_reps.rds"))

# Per-replicate diagnostics: theta alone and all four parameters
res <- res |>
  mutate(
    rhat_theta = map_dbl(diag, \(d) d$rhat[d$variable == "b_zmpd"]),
    ess_theta  = map_dbl(diag, \(d) min(d$ess_bulk[d$variable == "b_zmpd"],
                                       d$ess_tail[d$variable == "b_zmpd"])),
    rhat_max   = map_dbl(diag, \(d) max(d$rhat)),
    ess_min    = map_dbl(diag, \(d) min(d$ess_bulk, d$ess_tail)),
    clean_theta = rhat_theta < 1.01 & ess_theta > 400,
    clean_all   = rhat_max < 1.01 & ess_min > 400 & divergences == 0,
    excl = lo > 0 | hi < 0
  )

rate <- \(x) tibble(n = length(x), excl_zero = mean(x), mcse = sqrt(mean(x) * (1 - mean(x)) / length(x)))
summary_brms <- bind_rows(
  res |> group_by(design) |> reframe(subset = "all replicates", rate(excl)),
  res |> filter(clean_theta) |> group_by(design) |> reframe(subset = "theta converged", rate(excl)),
  res |> filter(clean_all) |> group_by(design) |> reframe(subset = "all parameters converged, no divergences", rate(excl))
) |> arrange(design)

diag_brms <- res |>
  select(design, rep, divergences, diag) |>
  unnest(diag) |>
  group_by(design, variable) |>
  summarise(rhat_max = max(rhat), rhat_gt_1.01 = sum(rhat > 1.01),
            ess_bulk_min = min(ess_bulk), ess_tail_min = min(ess_tail), .groups = "drop") |>
  left_join(res |> group_by(design) |> summarise(reps_with_divergences = sum(divergences > 0),
                                                 median_divergences = median(divergences)),
            by = "design")

# Do non-converged replicates differ in theta's interval?
by_clean <- res |> group_by(design, clean_all) |>
  summarise(n = n(), excl_zero = mean(excl), median_width = median(hi - lo), .groups = "drop")

print(summary_brms); print(diag_brms, n = Inf); print(by_clean)
write_csv(summary_brms, here("output", "tables", "design_sim_brms.csv"))
write_csv(diag_brms, here("output", "tables", "design_sim_brms_diagnostics.csv"))
