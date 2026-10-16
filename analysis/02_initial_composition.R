# Step 2: phylogenetic metrics fixed a priori (one per dimension of Tucker et
# al. 2017: PD = richness, MPD = divergence, VPD = regularity), computed for
# the initial (stocked, design-based) and final (realized, post-treatment)
# compositions.
#
# Initial-composition metrics are strongly correlated with the richness
# treatment (r(PD, div) = 0.99; output/tables/predictor_correlations.csv), so
# each metric is standardized within richness levels: slopes measure the effect
# of compositional differences at fixed richness.
#
# For each response we also fit a baseline model without phylogenetic
# predictors and report the proportion of among-composition variance
# (sd(comp)^2) explained by each metric.

library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)
library(here)

options(mc.cores = 4, brms.backend = "cmdstanr")
source(here("R", "brms_helpers.R"))

# ---- Data ----
bef <- read_csv(here("data", "derived", "bef_long.csv"), show_col_types = FALSE)

metrics <- c("PD", "MPD", "VPD")
comp_preds <- bef |>
  distinct(comp, div, across(c(paste0(metrics, "_init"), paste0(metrics, "_comp")))) |>
  rename_with(\(x) str_replace(x, "_comp$", "_final")) |>
  mutate(across(c(ends_with("_init"), ends_with("_final")), \(x) z_within(x, div), .names = "zw_{.col}"))

d <- bef |>
  left_join(select(comp_preds, comp, starts_with("zw_")), by = "comp") |>
  prep_data()

predictors <- c(paste0("zw_", metrics, "_init"), paste0("zw_", metrics, "_final"))

# Share of each metric's variance that lies within richness levels
within_share <- comp_preds |>
  pivot_longer(c(ends_with("_init"), ends_with("_final")) & !starts_with("zw_"),
               names_to = "metric") |>
  group_by(metric) |>
  summarise(r2_div = summary(lm(value ~ factor(div)))$r.squared, .groups = "drop") |>
  mutate(within_div_share = 1 - r2_div)
print(within_share)
write_csv(within_share, here("output", "tables", "step2_within_div_variance.csv"))

# ---- Fits ----
base <- map(responses, \(r) fit_uni(r, NULL, data = d, sigma_rhs = sigma98[[r]]))
fits <- expand_grid(resp = names(responses), x = predictors) |>
  mutate(fit = map2(resp, x, \(r, x) fit_uni(responses[[r]], x, data = d, sigma_rhs = sigma98[[responses[[r]]]])),
         name = paste("uni", resp, x, sep = "_"))

# ---- Diagnostics ----
diag <- bind_rows(
  imap_dfr(base, \(f, r) diagnose(f, paste0("uni_", r, "_base"))),
  pmap_dfr(list(fits$fit, fits$name), diagnose)
)
print(diag, n = Inf)
write_csv(diag, here("output", "tables", "step2_diagnostics.csv"))
walk2(fits$fit, fits$name, save_ppc)

# ---- Slopes ----
slopes <- pmap_dfr(list(fits$fit, fits$resp, fits$x), \(f, r, x)
  summarise_slopes(slope_draws(f, x), "univariate", r, x))
write_csv(slopes, here("output", "tables", "step2_slopes.csv"))

# Draws of the average within-richness slope, for Fig. 3
dir.create(here("output", "draws"), showWarnings = FALSE)
pmap_dfr(list(fits$fit, fits$resp, fits$x), \(f, r, x) {
  v <- slope_draws(f, x)$average
  tibble(response = r, predictor = x, .draw = seq_along(v), value = v)
}) |> saveRDS(here("output", "draws", "slope_draws_t98.rds"))

# ---- Among-composition variance explained ----
sd_comp <- \(f) as_draws_df(f)$sd_comp__Intercept
set.seed(SEED)
r2_comp <- pmap_dfr(list(fits$fit, fits$resp, fits$x), \(f, r, x) {
  s_base <- sd_comp(base[[r]])
  s_x <- sd_comp(f)
  v <- 1 - s_x^2 / sample(s_base^2, length(s_x), replace = TRUE)
  tibble(response = r, predictor = x,
         sd_comp_base = median(s_base), sd_comp_x = median(s_x),
         r2_comp_median = median(v),
         r2_comp_ratio_of_means = 1 - mean(s_x^2) / mean(s_base^2),
         q2.5 = quantile(v, 0.025), q97.5 = quantile(v, 0.975))
})
write_csv(r2_comp, here("output", "tables", "step2_r2_composition.csv"))

print(slopes |> filter(div == "average") |> mutate(across(where(is.numeric), \(x) round(x, 2))), n = Inf, width = Inf)
print(r2_comp |> mutate(across(where(is.numeric), \(x) round(x, 2))), n = Inf, width = Inf)
