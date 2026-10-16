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
# Outputs: output/tables/design_sim_brms.csv, design_sim_brms_diagnostics.csv

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
    tibble(
      design = nm, rep = r,
      est = median(dr$b_zmpd),
      lo = quantile(dr$b_zmpd, 0.025), hi = quantile(dr$b_zmpd, 0.975),
      max_rhat = max(rhat(fit), na.rm = TRUE),
      min_ess_bulk = min(summarise_draws(dr, ess_bulk)$ess_bulk, na.rm = TRUE),
      divergences = sum(np$Value[np$Parameter == "divergent__"])
    )
  })
}

res <- map_dfr(names(designs), run_design)

summary_brms <- res |>
  group_by(design) |>
  summarise(excl_zero = mean(lo > 0 | hi < 0),
            mcse = sqrt(excl_zero * (1 - excl_zero) / n()),
            median_width = median(hi - lo), .groups = "drop")
diag_brms <- res |>
  group_by(design) |>
  summarise(max_rhat = max(max_rhat), min_ess_bulk = min(min_ess_bulk),
            reps_with_divergences = sum(divergences > 0), .groups = "drop")

print(summary_brms); print(diag_brms)
write_csv(summary_brms, here("output", "tables", "design_sim_brms.csv"))
write_csv(diag_brms, here("output", "tables", "design_sim_brms_diagnostics.csv"))
