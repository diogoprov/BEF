# Sensitivity of the T99 primary models (within-cell MPD, initial composition)
# to the variance structure. Grouped posterior predictive checks show that the
# SD of productivity and respiration within some sampling periods is not
# reproduced even with sigma ~ period; this script checks whether alternative
# structures fix it and whether the focal slope depends on them.

library(tidyverse)
library(brms)
library(posterior)
library(here)

options(mc.cores = 4, brms.backend = "cmdstanr")
source(here("R", "brms_helpers.R"))
source(here("R", "t99_data.R"))

d <- prep_t99()
fixed <- "0 + cell + zw_MPD_init + ftime + dist + dist:post"
pri_sig <- c(priors_uni[priors_uni$class != "sigma", ], priors_sigma)

variants <- tribble(
  ~variant, ~random, ~sigma,
  "current: sigma ~ period, AR(1)", "(1 | comp) + (1 | tank) + ar(time = time, gr = tank, p = 1)", "ftime",
  "sigma ~ period x pulse, AR(1)", "(1 | comp) + (1 | tank) + ar(time = time, gr = tank, p = 1)", "ftime * dist",
  "sigma ~ period, no AR", "(1 | comp) + (1 | tank)", "ftime",
  "sigma ~ period, AR(1), + composition x period", "(1 | comp) + (1 | comp:ftime) + (1 | tank) + ar(time = time, gr = tank, p = 1)", "ftime"
)

check <- function(f, y) {
  set.seed(SEED)
  yr <- posterior_predict(f, ndraws = 300)
  yo <- f$data[[y]]
  g <- f$data$ftime
  p_sd <- map_dbl(levels(g), \(l) mean(apply(yr[, g == l], 1, sd) >= sd(yo[g == l])))
  p_m <- map_dbl(levels(g), \(l) mean(rowMeans(yr[, g == l]) >= mean(yo[g == l])))
  b <- as_draws_df(f)$b_zw_MPD_init
  tibble(sd_outside = sum(p_sd < 0.05 | p_sd > 0.95), mean_outside = sum(p_m < 0.05 | p_m > 0.95),
         periods = nlevels(g), slope = median(b), q2.5 = quantile(b, 0.025), q97.5 = quantile(b, 0.975))
}

res <- expand_grid(resp = c("prod", "resp"), variants) |>
  pmap_dfr(\(resp, variant, random, sigma) {
    y <- responses99[[resp]]
    pri <- if (grepl("ar(", random, fixed = TRUE)) pri_sig else pri_sig[pri_sig$class != "ar", ]
    f <- brm(bf(as.formula(paste(y, "~", fixed, "+", random)), as.formula(paste("sigma ~", sigma))),
             data = d, family = gaussian(), prior = pri, iter = ITER, chains = 4, seed = SEED,
             control = CONTROL, refresh = 0,
             file = here("output", "models", paste0("t99var_", resp, "_", make.names(variant))),
             file_refit = "on_change")
    bind_cols(tibble(response = resp, variant = variant), check(f, y), diagnose(f, variant) |> select(max_rhat, divergences))
  })

write_csv(res, here("output", "tables", "t99_variance_sensitivity.csv"))
print(res |> mutate(across(where(is.numeric), \(x) round(x, 3))), width = Inf)
