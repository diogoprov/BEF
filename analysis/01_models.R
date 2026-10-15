# Bayesian mixed models: phylogenetic diversity and ecosystem functioning
#
# Design: 84 tanks = 21 compositions x 4 tanks; 6 sampling periods per tank;
# richness treatment `div` (3, 9, 15 initial species) among compositions;
# acid pulse (`dist`) crossed with composition (2 tanks each), applied after
# sampling period 2. PD metrics are composition-level predictors.
#
# Model for each response (standardized; resp log-transformed; prod on the raw
# scale because log(prod) is left-skewed; dc with Student-t errors because of
# heavy tails, as indicated by posterior predictive checks):
#   y ~ 0 + div + div:zX + time + dist + dist:post
#       + (1 | comp) + (1 | tank) + AR(1) within tank
# `div:zX` gives one slope per richness level (taxonomic x phylogenetic
# diversity interaction); the average slope is derived from the draws.
#
# Multiple responses are handled in two ways:
#   (a) univariate models with 98.33% credible intervals (Bonferroni-type
#       alpha = 0.05/3, as suggested by A. Zuur & E. Ieno);
#   (b) a multivariate model (Student-t for all responses, since residual
#       correlations require a common family) with residual and
#       random-effect correlations.

library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)
library(here)

options(mc.cores = 4, brms.backend = "cmdstanr")
dir.create(here("output", "models"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("figures", "ppc"), showWarnings = FALSE, recursive = TRUE)

z <- \(x) as.numeric(scale(x))
SEED <- 20261015
ITER <- 2000
ITER_MV <- 4000
CONTROL <- list(adapt_delta = 0.95)

# ---- Data ----
bef <- read_csv(here("data", "derived", "bef_long.csv"), show_col_types = FALSE)

comp_preds <- bef |>
  distinct(comp, PD = PD_comp, MPD = MPD_comp, MNTD = MNTD_comp) |>
  mutate(zPD = z(PD), zMPD = z(MPD), zMNTD = z(MNTD))  # scaled across the 21 compositions

d <- bef |>
  left_join(select(comp_preds, comp, zPD, zMPD, zMNTD), by = "comp") |>
  mutate(
    y_prod = z(prod),
    y_resp = z(log(resp)),
    y_dc   = z(dc),
    div    = factor(div, levels = c(1, 9, 15), labels = c("S3", "S9", "S15")),
    ftime  = factor(time),
    comp   = factor(comp),
    tank   = factor(tank)
  )

responses <- c(prod = "y_prod", resp = "y_resp", dc = "y_dc")
families <- list(y_prod = gaussian(), y_resp = gaussian(), y_dc = student())
predictors <- c("zMPD", "zPD", "zMNTD")

# ---- Priors (standardized scale) ----
priors_uni <- c(
  prior(normal(0, 1), class = b),
  prior(exponential(2), class = sd),
  prior(exponential(1), class = sigma),
  prior(normal(0, 0.5), class = ar)
)
prior_nu <- prior(gamma(2, 0.1), class = nu)

rhs <- function(x, ar = TRUE) {
  paste0(
    "0 + div + div:", x, " + ftime + dist + dist:post + (1 | comp) + (1 | tank)",
    if (ar) " + ar(time = time, gr = tank, p = 1)" else ""
  )
}

fit_uni <- function(resp, x, ar = TRUE) {
  name <- paste("uni", resp, x, if (ar) "ar1" else "noar", sep = "_")
  brm(
    bf(paste(resp, "~", rhs(x, ar))),
    data = d, family = families[[resp]],
    prior = c(if (ar) priors_uni else priors_uni[priors_uni$class != "ar", ],
              if (families[[resp]]$family == "student") prior_nu),
    iter = ITER, chains = 4, seed = SEED, control = CONTROL,
    file = here("output", "models", name), file_refit = "on_change"
  )
}

# ---- Diagnostics helper ----
diagnose <- function(fit, name) {
  s <- summarise_draws(as_draws_df(fit), "rhat", "ess_bulk", "ess_tail")
  np <- nuts_params(fit)
  tibble(
    model = name,
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(s$ess_tail, na.rm = TRUE),
    divergences = sum(np$Value[np$Parameter == "divergent__"]),
    max_treedepth_hits = sum(np$Value[np$Parameter == "treedepth__"] >= 10)
  )
}

# Slopes of predictor x per richness level, plus their average
slope_draws <- function(fit, x, resp_prefix = "") {
  dr <- as_draws_df(fit)
  lv <- levels(d$div)
  cols <- paste0("b_", resp_prefix, "div", lv, ":", x)
  out <- map(cols, \(cl) dr[[cl]]) |> set_names(lv)
  out$average <- rowMeans(do.call(cbind, out))
  out
}

summarise_slopes <- function(sl, model, response, x) {
  imap_dfr(sl, \(v, k) tibble(
    model = model, response = response, predictor = x, div = k,
    mean = mean(v),
    q2.5 = quantile(v, 0.025), q97.5 = quantile(v, 0.975),
    q0.83 = quantile(v, 0.05 / 6), q99.17 = quantile(v, 1 - 0.05 / 6),
    p_gt0 = mean(v > 0)
  ))
}

save_ppc <- function(fit, name, resp = NULL) {
  p <- pp_check(fit, ndraws = 100, resp = resp) + ggtitle(name)
  ggsave(here("figures", "ppc", paste0(name, ".png")), p, width = 5, height = 3.5, dpi = 150)
}

# ---- Step 1: is AR(1) needed? (MPD models, with vs without AR) ----
ar_check <- imap_dfr(responses, \(r, nm) {
  f_ar <- fit_uni(r, "zMPD", ar = TRUE)
  f_no <- fit_uni(r, "zMPD", ar = FALSE)
  ar_dr <- as_draws_df(f_ar)$`ar[1]`
  lc <- loo_compare(list(ar1 = loo(f_ar), noar = loo(f_no)))
  tibble(
    response = nm,
    ar1_mean = mean(ar_dr), ar1_q2.5 = quantile(ar_dr, 0.025), ar1_q97.5 = quantile(ar_dr, 0.975),
    best_loo = rownames(lc)[1], elpd_diff = lc[2, "elpd_diff"], se_diff = lc[2, "se_diff"]
  )
})
write_csv(ar_check, here("output", "tables", "ar1_check.csv"))
print(ar_check)

# ---- Step 2: univariate models (AR1) for each response x predictor ----
uni <- expand_grid(resp = names(responses), x = predictors) |>
  mutate(fit = map2(resp, x, \(r, x) fit_uni(responses[[r]], x, ar = TRUE)),
         name = paste("uni", resp, x, sep = "_"))

diag_uni <- pmap_dfr(list(uni$fit, uni$name), diagnose)
slopes_uni <- pmap_dfr(list(uni$fit, uni$resp, uni$x), \(f, r, x)
  summarise_slopes(slope_draws(f, x), "univariate", r, x))
walk2(uni$fit, uni$name, save_ppc)

# Disturbance and richness effects from the MPD models (for reporting)
fixef_uni <- pmap_dfr(list(uni$fit, uni$resp, uni$x), \(f, r, x)
  as_tibble(fixef(f, probs = c(0.025, 0.975, 0.05 / 6, 1 - 0.05 / 6)), rownames = "term") |>
    mutate(response = r, predictor = x, .before = 1))

# ---- Step 3: multivariate MPD model ----
mv_bf <- function(resp) {
  bf(paste(resp, "~ 0 + div + div:zMPD + ftime + dist + dist:post",
           "+ (1 | c | comp) + (1 | t | tank) + ar(time = time, gr = tank, p = 1)"))
}
priors_mv <- c(
  prior(normal(0, 1), class = b, resp = yprod),
  prior(normal(0, 1), class = b, resp = yresp),
  prior(normal(0, 1), class = b, resp = ydc),
  prior(exponential(2), class = sd, resp = yprod),
  prior(exponential(2), class = sd, resp = yresp),
  prior(exponential(2), class = sd, resp = ydc),
  prior(exponential(1), class = sigma, resp = yprod),
  prior(exponential(1), class = sigma, resp = yresp),
  prior(exponential(1), class = sigma, resp = ydc),
  prior(normal(0, 0.5), class = ar, resp = yprod),
  prior(normal(0, 0.5), class = ar, resp = yresp),
  prior(normal(0, 0.5), class = ar, resp = ydc),
  prior(gamma(2, 0.1), class = nu),  # multivariate t: one shared nu
  prior(lkj(2), class = cor),
  prior(lkj(2), class = rescor)
)
fit_mv <- brm(
  mv_bf("y_prod") + mv_bf("y_resp") + mv_bf("y_dc") + set_rescor(TRUE),
  data = d, family = student(), prior = priors_mv,  # rescor requires one family
  iter = ITER_MV, chains = 4, seed = SEED, control = CONTROL,
  file = here("output", "models", "mv_zMPD"), file_refit = "on_change"
)

diag_mv <- diagnose(fit_mv, "mv_zMPD")
slopes_mv <- imap_dfr(c(prod = "yprod", resp = "yresp", dc = "ydc"), \(rp, nm)
  summarise_slopes(slope_draws(fit_mv, "zMPD", paste0(rp, "_")), "multivariate", nm, "zMPD"))
iwalk(c(prod = "yprod", resp = "yresp", dc = "ydc"), \(rp, nm) save_ppc(fit_mv, paste0("mv_zMPD_", nm), resp = rp))

fixef_mv <- as_tibble(fixef(fit_mv, probs = c(0.025, 0.975)), rownames = "term") |>
  mutate(model = "mv_zMPD", .before = 1)

# ---- Write outputs ----
bind_rows(diag_uni, diag_mv) |> write_csv(here("output", "tables", "diagnostics.csv"))
bind_rows(slopes_uni, slopes_mv) |> write_csv(here("output", "tables", "pd_slopes.csv"))
write_csv(fixef_uni, here("output", "tables", "fixef_univariate.csv"))
write_csv(fixef_mv, here("output", "tables", "fixef_multivariate.csv"))
print(summary(fit_mv)$rescor_pars)

print(bind_rows(diag_uni, diag_mv), n = Inf)
print(bind_rows(slopes_uni, slopes_mv) |> mutate(across(where(is.numeric), \(x) round(x, 3))), n = Inf, width = Inf)
