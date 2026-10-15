# Prior predictive checks and prior sensitivity (BARG; Kruschke 2021) for the
# primary models: within-richness/within-cell MPD (initial composition) in T98
# (analysis/02_initial_composition.R) and T99 (analysis/11_models_t99.R, part a).
#
# (1) Prior predictive: each model refitted with sample_prior = "only"; the
#     prior predictive distribution of the standardized response is compared
#     with the observed one.
# (2) Prior sensitivity: power-scaling sensitivity (priorsense) of the focal
#     parameters in the fitted models.

library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)
library(priorsense)
library(here)

options(mc.cores = 4, brms.backend = "cmdstanr")
source(here("R", "brms_helpers.R"))
source(here("R", "t99_data.R"))
dir.create(here("figures", "prior_predictive"), showWarnings = FALSE, recursive = TRUE)

# ---- Data ----
bef <- read_csv(here("data", "derived", "bef_long.csv"), show_col_types = FALSE)
d98 <- bef |>
  left_join(bef |> distinct(comp, div, MPD_init) |>
              mutate(zw_MPD_init = z_within(MPD_init, div)) |> select(comp, zw_MPD_init),
            by = "comp") |>
  prep_data()
d99 <- prep_t99()

rhs_a99 <- \(resp) paste0("0 + cell + zw_MPD_init + ftime + dist + dist:post", random_for(resp))

models <- bind_rows(
  tibble(exp = "T98", resp = names(responses), y = unname(responses),
         name = paste0("uni_", unname(responses), "_zw_MPD_init_ar1", c("_sig", "_sig", ""))),
  tibble(exp = "T99", resp = names(responses99), y = unname(responses99),
         name = paste0("t99a_", names(responses99), "_zw_MPD_init_sig"))
)

get_fit <- function(exp, resp, y, name, prior_only = FALSE) {
  nm <- if (prior_only) paste0(name, "_prioronly") else name
  if (exp == "T98") {
    f <- fit_uni(y, "zw_MPD_init", data = d98, sigma_rhs = sigma98[[y]])
  } else {
    f <- fit_uni(y, data = d99, name = name, rhs_override = rhs_a99(resp), sigma_rhs = sigma99[[resp]])
  }
  if (!prior_only) return(f)
  update(f, sample_prior = "only", recompile = TRUE, seed = SEED, refresh = 0,
         file = here("output", "models", nm), file_refit = "on_change")
}

# ---- (1) Prior predictive ----
prior_pred <- pmap_dfr(models, \(exp, resp, y, name) {
  fp <- get_fit(exp, resp, y, name, prior_only = TRUE)
  yrep <- posterior_predict(fp, ndraws = 200)
  yobs <- fp$data[[y]]
  p <- ppc_dens_overlay(yobs, yrep[1:100, ]) + coord_cartesian(xlim = c(-10, 10)) +
    ggtitle(paste(exp, resp, "prior predictive"))
  ggsave(here("figures", "prior_predictive", paste0(exp, "_", resp, ".png")), p,
         width = 5, height = 3.5, dpi = 150)
  tibble(exp = exp, response = resp,
         obs_sd = sd(yobs), prior_pred_sd_median = median(apply(yrep, 1, sd)),
         prior_pred_q99_abs = quantile(abs(yrep), 0.99),
         share_within_obs_range = mean(yrep >= min(yobs) & yrep <= max(yobs)))
})
write_csv(prior_pred, here("output", "tables", "prior_predictive_summary.csv"))
print(prior_pred)

# ---- (2) Prior sensitivity (power scaling) ----
sens <- pmap_dfr(models, \(exp, resp, y, name) {
  f <- get_fit(exp, resp, y, name)
  vars <- variables(f)
  focal <- vars[str_detect(vars, "zw_MPD_init|^b_dist|^sd_|^sigma$|^b_sigma_Intercept|^ar\\[")]
  s <- powerscale_sensitivity(f, variable = focal)
  as_tibble(s) |> mutate(exp = exp, response = resp, .before = 1)
})
write_csv(sens, here("output", "tables", "prior_sensitivity.csv"))
print(sens, n = Inf, width = Inf)
